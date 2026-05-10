"""Base classes for fitness models."""
from typing import Any
from collections.abc import Callable, Iterable, Mapping
from abc import ABC, abstractmethod

from carabiner import print_err
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
        weights: ArrayLike | None = None, 
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
        valid: ArrayLike | None = None, 
        weights: ArrayLike | None = None, 
        param_names: str | None = None, 
        **kwargs
    ) -> dict[str, float | int]:
        n_orig = len(y)
        if valid is None:
            valid = np.ones_like(y, dtype=bool)
        y = y[valid]
        x = x[valid]
        if weights is not None:
            weights = np.asarray(weights)[valid]
        weights = self.calculate_weights(y, weights)
        betas, cis, ses, ps, preds, other = self._fit(y, x, weights=weights, param_names=param_names, **kwargs)

        # pad preds and y back to original shape for consistent stacking
        preds_full = np.full(n_orig, np.nan)
        preds_full[valid] = preds
        x_out = np.full(n_orig, np.nan)
        x_out[valid] = x
        y_out = np.full(n_orig, np.nan)
        y_out[valid] = y

        result = betas | {
            f"{k}_p": v 
            for k, v in ps.items()
        } | {
            f"{k}_se": v for k, v in ses.items()
        } | {
            f"{k}_ci_{j}": _v 
            for k, v in cis.items() 
            for j, _v in zip(["low", "high"], v)
        } | {
            "nobs": y.shape[0],
        } | other
        return self._fitness_transform(result), (x_out, y_out, preds_full)

    def fit(
        self, 
        Y: ArrayLike, 
        x: ArrayLike, 
        valid: ArrayLike | None = None, 
        weights: ArrayLike | None = None, 
        param_names: str | None = None,
        groups: ArrayLike | None = None,
        min_obs: int = 3,
        weight_kwargs: Mapping[str, Any] | None = None,
        **kwargs
    ) -> list[dict[str, float | int]]:
            
        if valid is None:
            valid = np.ones(Y.shape[0], dtype=bool)
        finite_X = np.isfinite(x)
        finite_Y = np.isfinite(Y)

        if groups is not None:
            unique_groups = sorted(set(groups))

        results = []
        base_result = {"fit_status": "fail"}
        if weight_kwargs is None:
            weight_kwargs = {}
        if weights is None:
            if len(weight_kwargs) > 0:
                print_err(f"[INFO] Calculating observation weights: {weight_kwargs}")
                weights = self.calculate_weights_matrix(Y, weights=None, **weight_kwargs)
            else:
                print_err("[INFO] Not calculating observation weights")
                weights = [None] * Y.shape[0]
        else:
            print_err(f"[INFO] Calculating observation weights: {weights=}, {weight_kwargs}")
            weights = self.calculate_weights_matrix(Y, weights, **weight_kwargs)
        if Y.shape[0] > 10:
            from tqdm.auto import tqdm
            Y = tqdm(Y, desc="Fitting models")
        _preds = []
        for i, (y, finite_y, w, _valid) in enumerate(zip(
            Y,
            finite_Y,
            weights,
            valid,
        )):
            this_base_result = base_result | {"i": i}
            this_valid = _valid & finite_y & finite_X
            if w is not None:
                this_valid = this_valid & np.isfinite(w) & (w > 0.)

            if this_valid.sum() < min_obs:
                results.append(this_base_result)
                _preds.append((np.array([]), np.full(x.shape, np.nan), np.full(x.shape, np.nan)))
                continue
            if groups is None:
                result, preds = self.fit_obs(
                    y, 
                    x, 
                    valid=this_valid, 
                    weights=w, 
                    param_names=param_names,
                    **kwargs,
                )
            else:
                result = {}
                preds = (np.empty_like(x), np.empty_like(y), np.empty_like(y))
                for _group in unique_groups:
                    group_mask = groups == _group
                    _x, _y = x[group_mask], y[group_mask]
                    _result, _preds_new = self.fit_obs(
                        _y,
                        _x,
                        valid=this_valid[group_mask],
                        weights=None if w is None else w[group_mask], 
                        param_names=param_names,
                        **kwargs,
                    )
                    result |= {f"{key}@{_group}": v for key, v in _result.items()}
                    for _p, new_pred in zip(preds, _preds_new):
                        _p[group_mask] = new_pred
            results.append(this_base_result | result | {"fit_status": "ok"})
            _preds.append(preds)
        _preds = (
            np.concatenate([_x for _x, *_ in _preds]),
            np.stack([_y for _, _y, _ in _preds]),
            np.stack([_p for _, _, _p in _preds]),
        )
        return results, _preds


class LinearModel(Model):

    """Linear models using statsmodels."""

    def __init__(self, *args, method: Callable, **kwargs):
        super().__init__(*args, **kwargs)
        self.method = method
    
    def _fit(
        self, 
        y: ArrayLike, 
        x: ArrayLike,
        param_names: Iterable[str] | None = None, 
        method: Callable | None = None, 
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
        preds = res.predict(x)
        return betas, cis, ses, ps, preds, {"rsq": res.rsquared}

    
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
        param_names: Iterable[str] | None = None, 
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
                # method="dogbox",
            )
        except RuntimeError:
            betas = np.full(len(param_names), np.nan)
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
        preds = model_fn(x, *(v for _, v in betas.items()))
        return betas, cis, ses, ps, preds, {"dof": dof} | {f"{p}_init": v for p, v in zip(param_names, init_params)}
