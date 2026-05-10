
from typing import TYPE_CHECKING, Any, Iterable, Optional, Union
from abc import ABC, abstractmethod

if TYPE_CHECKING:
    from anndata import AnnData
else:
    AnnData = Any
from numpy.typing import ArrayLike

from .base import Model
from .linear import OLSModel, WLSModel
from .nonlinear import HillFitnessModel
from ..transforms import compute_log_ratios


class AnnDataModel(Model, ABC):

    def fit(
        self,
        adata: AnnData,
        Y_key: str = "__log_ratio__",
        X_key: str = "__log_expansion__",
        valid_keys: Optional[Union[str, Iterable[str]]] = None,
        pseudocount: float = 0.,
        volume_column: Optional[str] = None,
        _index_col: str = "__index__",
        **kwargs
    ) -> AnnData:
        """Fit fitness model per strain.

        Parameters
        ==========
        adata       : must have layers['log_ratio'] and var['__expansion__']
        model       : Model instance, defaults to self
        dispersion  : (n_strains,) pre-estimated dispersions; estimated from 
                    replicates if None
        min_obs     : minimum valid observations to attempt fit
        pseudocount : added to raw counts before weight calculation

        """
        import pandas as pd
        import numpy as np

        if hasattr(self, "_name"):
            name = self._name or "AnnDataModel"
        else:
            name = "AnnDataModel"
        if Y_key not in adata.layers:
            adata = compute_log_ratios(
                adata,
                pseudocount=pseudocount,
                volume_column=volume_column,
            )

        if valid_keys is not None:
            if isinstance(valid_keys, str):
                valid_keys = [valid_keys]
            valid = np.stack([adata.obs[k].values for k in valid_keys], axis=-1)
            valid = np.all(valid, axis=1)
        else:
            valid = np.ones((adata.X.shape[0],), dtype=bool)

        results, (x, y, preds) = super().fit(
            Y=adata.layers[Y_key], 
            x=adata.var[X_key].values, 
            valid=valid,
            **kwargs
        )
        adata.layers[f"{name}:predicted"] = preds
        adata.layers[f"{name}:residual"] = y - preds

        results_with_index = []
        for i, (idx, res) in enumerate(zip(
            adata.obs.index, 
            results,
        )):
            results_with_index.append(res | {_index_col: idx})

        results_df = pd.DataFrame(results_with_index, index=adata.obs.index)
        results_df = results_df.rename(columns={
            col: f"{name}:{col}" 
            for col in results_df 
            if col != _index_col
        })
        # transparent writeback — any keys model returns land in obs
        adata.obs = adata.obs.merge(
            results_df.set_index(_index_col), 
            how="left", 
            validate="one_to_one",
            left_index=True, 
            right_index=True,
        )
        adata.uns["models_fitted"].append(name)
        return adata


class AnnDataWLSModel(AnnDataModel, WLSModel):
    def fit(
        self, 
        adata, 
        *args, 
        weight_kwargs: None = None, 
        concentration_key: str = "__inducer__", 
        **kwargs
    ):
        return super().fit(
            adata,
            *args,
            groups=(
                adata.var[concentration_key].values 
                if adata.uns["concentration_column"] and adata.var[concentration_key].nunique() > 1
                else None
            ),
            weight_kwargs={
                "raw": adata.X,
                "control_mask": adata.obs["__is_reference__"].values,
                "groups": adata.var["__culture_index__"].values,
            },
            **kwargs,
        )


class AnnDataOLSModel(AnnDataModel, OLSModel):
    def fit(
        self, 
        adata, 
        *args, 
        weights: None = None, 
        concentration_key: str = "__inducer__", 
        **kwargs
    ):
        return super().fit(
            adata,
            *args,
            groups=adata.var[concentration_key].values if adata.uns["concentration_column"] else None,
            weights=None,
            weight_kwargs={},
            **kwargs,
        )


class AnnDataHillModel(AnnDataModel, HillFitnessModel):
    def fit(
        self, 
        adata, 
        *args, 
        concentration: str, 
        weight_kwargs: None = None, 
        concentration_key: str = "__inducer__", 
        **kwargs
    ):
        return super().fit(
            adata,
            *args,
            concentration=adata.var[concentration_key].values,
            weight_kwargs={
                "raw": adata.X,
                "control_mask": adata.obs["__is_reference__"].values,
                "groups": adata.var["__culture_index__"].values,
            },
            **kwargs,
        )
