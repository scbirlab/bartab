"""Non-linear models."""
from typing import Any, Callable, Dict, Iterable, Mapping, Optional
from functools import partial

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import expit as sigmoid

from .base import NonLinear

def hill_function(log_c, log_ic50, bottom=0., h=1., top=1.):
    return bottom + (top - bottom) * sigmoid(h * (log_ic50 - log_c))


class HillFitnessModel(NonLinear):
    """2-parameter Hill / log-logistic with bottom fixed at 0.
    
    Parameters: log(EC50) [location] and H [slope], fit on log concentration.
    Fitness = 1 - σ(H · log(c/EC50))
    """

    _name: str = "HillFitnessModel"

    @staticmethod
    def _model_fn(
        X: Iterable[ArrayLike],
        *args
    ):
        x, log_c = X
        return x * (1. - hill_function(log_c, *args))

    @staticmethod
    def _fitness_transform(results: Mapping[str, float]) -> Dict[str, float]:
        fn = np.exp
        param = "log_ic50"
        return results | {
            "ic50": fn(results[param]), 
            "ic50_high": fn(results[f"{param}_ci_high"]),
            "ic50_low": fn(results[f"{param}_ci_low"]),
        }

    @staticmethod
    def _per_concentration_fitness(y, x, weights, log_c) -> Dict[float, float]:
        import statsmodels.api as sm
        fitness_by_conc = {}
        for c_val in np.unique(log_c):
            mask = log_c == c_val
            if mask.sum() < 2:
                continue
            res = sm.WLS(
                y[mask], 
                x[mask], 
                weights=weights[mask],
            ).fit()
            fitness_by_conc[c_val] = 1. - res.params[0]
        return fitness_by_conc

    @staticmethod
    def _init_params(y, x, weights, log_c, fitness_by_conc):
        if len(fitness_by_conc) < 2:
            return np.median(log_c), 1.

        log_c_observed = sorted(fitness_by_conc)
        inhibition = 1. - np.array([
            fitness_by_conc[c] for c in log_c_observed
        ])
        # log_ec50: interpolate where inhibition crosses 0.5
        if (inhibition > .5).any() and (inhibition < .5).any():
            log_ec50_init = np.interp(.5, inhibition, log_c_observed)
        elif inhibition.mean() < .5:
            log_ec50_init = log_c_observed[-1]   # effect beyond measured range
        else:
            log_ec50_init = log_c_observed[0]    # effect below measured range

        # H: gradient near EC50, using logistic identity slope = H/4 at inflection
        gradients = np.diff(inhibition) / np.diff(log_c_observed)
        idx = np.argmin(np.abs(inhibition[:-1] - .5))
        slope = 4. * gradients[idx]
        h_init = np.clip(slope, .1, 10.)
        return log_ec50_init, h_init

    def _fit(
        self, 
        y: ArrayLike, 
        x: ArrayLike,
        concentration: ArrayLike,
        weights: ArrayLike,
        init_params: Optional[ArrayLike] = None, 
        param_names: Optional[Iterable[str]] = None,
        model_fn: Optional[Callable] = None,
        model_kwargs: Optional[Mapping[str, Any]] = None,
        n_params: int = 2,
        **kwargs
    ):
        if param_names is None:
            param_names = ["log_ic50", "bottom", "h", "top"][:n_params]
        zero_mask = concentration == 0.
        if zero_mask.any():
            valid = ~zero_mask
            y, x, weights, concentration = (
                arr[valid] for arr in (y, x, weights, concentration)
            )
        log_c = np.log(concentration)
        if init_params is None:
            fitness_by_conc = self._per_concentration_fitness(y, x, weights, log_c)
            init_params = self._init_params(y, x, weights, log_c, fitness_by_conc)[:n_params]
        if model_fn is None:
            model_fn = self._model_fn
        if model_kwargs is None:
            model_kwargs = {}
        return super()._fit(
            y=y, 
            x=(x, log_c), 
            model_fn=partial(model_fn, **model_kwargs),
            init_params=init_params, 
            weights=weights, 
            param_names=param_names,
            **kwargs
        )
