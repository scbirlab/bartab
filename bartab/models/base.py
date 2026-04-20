"""Base classes for fitness models."""
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Union
from abc import ABC, abstractmethod

import pandas as pd
import numpy as np
from numpy.typing import ArrayLike


class Model(ABC):

    _name: str = None

    def __init__(
        self,
        ci: float = .95
    ):
        self.ci_level = ci

    @abstractmethod
    def _fit(self, y: ArrayLike, x: ArrayLike, weights: ArrayLike, param_names: Iterable[str], **kwargs):
        ...

    @staticmethod
    def _fitness_transform(results: Mapping[str, float], **kwargs):
        return results

    @staticmethod
    def _calculate_weights(y: ArrayLike, **kwargs) -> np.ndarray:
        return np.ones_like(y)

    def calculate_weights(self, y: ArrayLike, weights=None, **kwargs) -> np.ndarray:
        if weights is None:
            return self._calculate_weights(y)
        else:
            # pass through
            return np.asarray(weights)

    def _calculate_weights_matrix(
        self, 
        Y: ArrayLike,
        **kwargs
    ) -> np.ndarray:
        return np.stack([
            self._calculate_weights(y, **kwargs) 
            for y in Y
        ], axis=0)

    def calculate_weights_matrix(
        self, 
        Y: ArrayLike, 
        weights: Optional[ArrayLike] = None, 
        **kwargs
    ) -> np.ndarray:
        if weights is None:
            return self._calculate_weights_matrix(Y, **kwargs)
        else:
            # pass through
            return np.asarray(weights)

    def fit_obs(
        self, 
        y: ArrayLike, 
        x: ArrayLike, 
        valid: Optional[ArrayLike] = None, 
        weights: Optional[ArrayLike] = None, 
        param_names: Optional[str] = None, 
        **kwargs
    ) -> Dict[str, Union[float, int]]:
        if valid is None:
            valid = np.ones_like(y)
        y = y[valid]
        x = x[valid]
        weights = self.calculate_weights(y, weights)
        betas, cis, ses, ps, other = self._fit(y, x, weights=weights, param_names=param_names, **kwargs)
        result = betas | {
            f"{k}_p": v 
            for k, v in ps.items()
        } | {
            f"{k}_se": v for k, v in ses.items()
        } | {
            f"{k}_ci_{j}": _v 
            for k, v in cis.items() 
            for j, _v in zip(["low", "high"], v)
        } | {"nobs": y.shape[0]} | other
        return self._fitness_transform(result)

    def fit(
        self, 
        Y: ArrayLike, 
        x: ArrayLike, 
        valid: Optional[ArrayLike] = None, 
        weights: Optional[ArrayLike] = None, 
        param_names: Optional[str] = None, 
        min_obs: int = 3,
        weight_kwargs: Optional[Mapping[str, ...]] = None,
        **kwargs
    ) -> List[Dict[str, Union[float, int]]]:
        if valid is None:
            valid = np.ones(Y.shape[0], dtype=bool)
        finite_X = np.isfinite(x)
        finite_Y = np.isfinite(Y)
        results = []
        base_result = {"fit_status": "fail"}
        if weight_kwargs is None:
            weight_kwargs = {}
        if weights is None:
            weights = [None] * Y.shape[0]
        else:
            weights = self.calculate_weights_matrix(Y, weights, **weight_kwargs)
        for i, (y, finite_y, w, _valid) in enumerate(zip(
            Y,
            finite_Y,
            weights,
            valid
        )):
            this_base_result = base_result | {"__i__": i}
            this_valid = _valid & finite_y & finite_X
            if w is not None:
                this_valid = this_valid & np.isfinite(w) & (w > 0.)

            if this_valid.sum() < min_obs:
                results.append(this_base_result)
                continue

            result = self.fit_obs(
                y, 
                x, 
                valid=this_valid, 
                weights=w, 
                param_names=param_names,
                **kwargs,
            )
            results.append(this_base_result | result | {"fit_status": "ok"})

        return results


class LinearModel(Model):

    """Linear models using statsmodels."""

    def __init__(self, *args, method: Callable, **kwargs):
        super().__init__(*args, **kwargs)
        self.method = method
    
    def _fit(
        self, 
        y: ArrayLike, 
        x: ArrayLike,
        param_names: Optional[Iterable[str]] = None, 
        method: Optional[Callable] = None, 
        **kwargs
    ):
        if method is None:
            method = self.method
        res = self.method(y, x, **kwargs).fit()
        if param_names is None:
            param_names = list(map(str, range(len(res.params))))
        betas = dict(zip(param_names, res.params))
        ses = dict(zip(param_names, res.bse))
        cis = dict(zip(param_names, res.conf_int(alpha=1. - self.ci_level)))
        ps = dict(zip(param_names, res.pvalues))
        return betas, cis, ses, ps, {"rsq": res.rsquared}

    
class NonLinear(Model):

    """Non-linear models using scipy.curve_fit.
    """

    _name: str = "NonLinear"
    
    def _fit(
        self, 
        y: ArrayLike, 
        x: ArrayLike, 
        model_fn: Callable,
        init_params: ArrayLike, 
        weights: ArrayLike, 
        param_names: Optional[Iterable[str]] = None, 
        **kwargs
    ):
        from scipy.optimize import curve_fit
        from scipy.stats import t as tdist
        if param_names is None:
            param_names = list(map(str, range(len(init_params))))
        dof = max(0, y.shape[0] - len(param_names))
        try:
            betas, pcov = curve_fit(
                model_fn,
                x,
                y,
                p0=init_params,
                sigma=1. / np.sqrt(weights),
                absolute_sigma=True,
            )
        except RuntimeError:
            betas = {
                k: np.nan for k in param_names
            }
            pcov = np.full(
                (len(param_names), len(param_names)), 
                np.nan,
            )

        betas = dict(zip(param_names, betas))
        ses = dict(zip(param_names, np.sqrt(np.diag(pcov))))
        cis = {
            k: tdist.interval(self.ci_level, dof, loc=b, scale=ses[k]) 
            for k, b in betas.items()
        }
        ps = {
            k: 2. * tdist.cdf(-abs(b / ses[k]), dof)
            for k, b in betas.items()
        }
        return betas, cis, ses, ps, {"dof": dof}
