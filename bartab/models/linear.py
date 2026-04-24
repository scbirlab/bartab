"""Fitness estimation."""
from typing import Callable, Dict, Iterable, Mapping, Optional

import numpy as np
from numpy.typing import ArrayLike

from .base import LinearModel


def _delta_method_weights(
    raw: ArrayLike,           # (n_strains, n_samples) raw counts + pseudocount
    control_mask: ArrayLike,    # (n_strains,)
    dispersion: ArrayLike,  # (n_strains,) per-strain alpha
) -> np.ndarray:
    """
    
    Var(log c_i - log c_wt) ≈ 1/c_i + α_i + 1/c_wt + α_wt → weights = 1/Var.
    
    Examples
    ========
    >>> import numpy as np; np.set_printoptions(legacy='1.25')
    >>> from bartab.models.linear import _delta_method_weights
    >>> raw = np.array([[500., 1_000., 2_000.],   # reference (wt)
    ...                 [100.,   150.,   200.]])   # mutant
    >>> control_mask = np.array([True, False])
    >>> dispersion = np.array([0.01, 0.05])
    >>> weights = _delta_method_weights(raw, control_mask, dispersion)
    >>> weights.shape  # (n_strains, n_samples)
    (2, 3)
    >>> (weights > 0).all()
    True
    
    """
    ref_counts = raw[control_mask, :].squeeze(axis=0)   # (n_samples,)
    ref_disp = dispersion[control_mask].squeeze()    # scalar
    var_y = (
        1. / raw + dispersion[:, None]              # (n_strains, n_samples)
        + 1. / ref_counts + ref_disp
    )
    return 1. / var_y


def _estimate_dispersion_mom(
    raw: ArrayLike, 
    groups: Optional[ArrayLike] = None
) -> np.ndarray:
    """Method-of-moments per-strain dispersion from within-culture replicate variance.
    
    Fallback when PyDESeq2 is unavailable or pool is too small.
    α̂ = (σ²/μ² - 1/μ)  clipped to [0, ∞)

    Examples
    ========
    >>> import numpy as np; np.set_printoptions(legacy='1.25')
    >>> from bartab.models.linear import _estimate_dispersion_mom
    >>> raw = np.array([[100., 200., 100.],   # high variance
    ...                 [100., 100., 100.]])  # zero variance
    >>> disp = _estimate_dispersion_mom(raw)
    >>> disp.shape
    (2,)
    >>> disp[1]  # constant counts → no extra-Poisson dispersion
    0.0
    >>> disp[0] > 0  # noisy counts → positive dispersion
    True

    """
    dispersions = []
    if groups is None:
        groups = np.ones((raw.shape[1]))
    for i in range(raw.shape[0]):
        means, _vars = [], []
        for group in np.unique(groups):
            vals = raw[i, groups == group]
            means.append(vals.mean())
            _vars.append(vals.var())
        mu = np.mean(means)
        dispersions.append(
            max(0., np.mean(_vars) / np.square(mu) - 1. / mu)
        )
    return np.asarray(dispersions)


class WLSModel(LinearModel):

    _name: str = "WLS"

    def __init__(self, *args, **kwargs):
        import statsmodels.api as sm
        super().__init__(*args, method=sm.WLS, **kwargs)
    
    def _calculate_weights_matrix(
        self,
        Y: ArrayLike,  # .fit() expects Y first argument
        raw: ArrayLike, 
        control_mask: ArrayLike,
        pseudocount: float = .5,
        groups: Optional[ArrayLike] = None
    ) -> np.ndarray:
        raw = raw + pseudocount
        dispersion = _estimate_dispersion_mom(
            raw,
            groups,
        )
        return _delta_method_weights(raw, control_mask, dispersion)  # (n_strains, n_samples)

    @staticmethod
    def _fitness_transform(results: Mapping[str, float]) -> Dict[str, float]:
        fn = lambda x: 1. - x
        param = "slope"
        return results | {
            "fitness": fn(results[param]), 
            "fitness_high": fn(results[f"{param}_ci_low"]),
            "fitness_low": fn(results[f"{param}_ci_high"]),
        }

    def _fit(self, 
        *args,
        param_names: Optional[Iterable[str]] = None,
        method: Optional[Callable] = None, 
        weights: Optional[ArrayLike] = None,
        **kwargs
    ):
        import statsmodels.api as sm
        if param_names is None:
            param_names = ["slope"]
        if weights is None:
            method = sm.OLS
            return super()._fit(
                *args, 
                param_names=param_names, 
                method=method, 
                **kwargs,
            )
        else:
            if method is None:
                method = self.method
            return super()._fit(
                *args, 
                param_names=param_names, 
                method=method,
                weights=weights,
                **kwargs,
            )


class OLSModel(WLSModel):
    """Unweighted baseline for benchmarking.
    
    Examples
    ========
    >>> import numpy as np; np.set_printoptions(legacy='1.25')
    >>> from bartab.models.linear import OLSModel
    >>> model = OLSModel()
    >>> Y = np.array([[0.0, -0.5, -1.0],   # declining strain
    ...               [0.0,  0.0,  0.0]])  # neutral strain
    >>> x = np.array([0., 1., 2.])
    >>> results, _ = model.fit(Y, x)
    >>> len(results)
    2
    >>> results[0]["fitness"]  # slope=-0.5 → fitness = 1-(-0.5) = 1.5
    1.5
    >>> results[1]["fitness"]  # slope=0 → fitness = 1.0
    1.0
    
    """

    _name: str = "OLS"

    def _fit(self, *args, weights=None, **kwargs):
        return super()._fit(*args, weights=None, **kwargs)
