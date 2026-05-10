"""Transforms on AnnData."""
from typing import TYPE_CHECKING, Any, Callable, Iterable, Literal, Optional, Union

if TYPE_CHECKING:
    from anndata import AnnData
else:
    AnnData = Any

from carabiner import print_err
import numpy as np


def _anndata_groupby(
    adata: AnnData,
    grouping: Union[str, Iterable[str]],
    axis: str = "var",
    **kwargs
):
    if isinstance(axis, int):
        axis = ("obs", "var")[axis]
    df = getattr(adata, axis)
    for group, group_df in df.groupby(grouping, **kwargs):
        idx = df.index.isin(group_df.index)
        yield group, group_df, idx


def compute_log_ratios(
    adata: AnnData,
    pseudocount: float = 0.,
    volume_column: Optional[str] = None,
    growth_column: Optional[str] = None,
    growth_type: Literal["density", "generations"] = "density",
    use_spike: bool = False,
    key: str = "__log_ratio__"
) -> AnnData:
    """Compute log-ratio layers for fitness regression.

    Adds layers:
        log_ratio  : log(c_i(t)/c_wt(t) * c_wt(0)/c_i(0))  [y-axis]
        expansion  : log(c_spike(t)/c_wt(t) * c_wt(0)/c_spike(0))  [x-axis]
                     with optional volume correction for adaptive sampling.
    
    """
    X = adata.X + pseudocount  # (n_strains, n_samples)

    t0_mask  = adata.var["__is_t0__"].values         # (n_samples,)
    ref_mask = adata.obs["__is_reference__"].values   # (n_strains,)
    spike_mask = adata.obs["__is_spike__"].values     # (n_strains,)

    n_ref = ref_mask.sum()
    if n_ref == 0:
        raise ValueError("No reference barcodes found")
    if n_ref > 0:
        print_err(f"Found {n_ref} reference barcodes: {', '.join(adata.obs.index[ref_mask].astype(str))}")
    if t0_mask.sum() == 0:
        raise ValueError("No t0 samples found")

    log_X = np.log(X)                                # (n_strains, n_samples)
    log_ref = np.log(np.sum(X[ref_mask, :], axis=0, keepdims=True))    # (1, n_samples,)

    # log(c_i / c_wt) at every sample
    log_ratio_to_ref = log_X - log_ref      # (n_strains, n_samples)
    log_ratio = log_ratio_to_ref.copy()

    # subtract per-culture t0 baseline
    for culture, _, c_idx in _anndata_groupby(adata, "__culture_index__"):
        t0_idx = c_idx & t0_mask
        if t0_idx.sum() == 0:
            raise ValueError(f"No t0 sample found for culture '{culture}'")
        log_ratio[:, c_idx] -= log_ratio_to_ref[:, t0_idx].mean(axis=1, keepdims=True)

    adata.layers[key] = log_ratio

    # expansion axis from spike
    if use_spike and spike_mask.any():
        n_spikes = spike_mask.sum()
        if n_spikes != 1:
            raise ValueError(f"Expected exactly 1 spike strain, found {n_spikes}")
        print_err(f"[INFO] Calculating culture expansion using spike-in")
        expansion = log_ratio.copy()[spike_mask, :]

        # volume correction for adaptive-volume sampling
        if volume_column is not None:
            if volume_column not in adata.var.columns:
                raise ValueError(f"Volume_column '{volume_column}' not in sample metadata")
            vols = adata.var[volume_column].values.astype(float)
            for culture, _, c_idx in _anndata_groupby(adata, "__culture_index__"):
                t0_idx = c_idx & t0_mask
                vol_t0_mean = vols[t0_idx].mean()
                expansion[:, c_idx] += (np.log(vol_t0_mean) - np.log(vols[c_idx]))[None, :]

    elif growth_column is not None:
        if growth_column not in adata.var.columns:
            raise ValueError(f"Column '{growth_column=}' not in sample metadata.")
        print_err(f"[INFO] Calculating culture expansion using {growth_column=}")
        growth = adata.var[growth_column].values
        if (growth <= 0).any():
            raise ValueError(
                f"'{growth_column=}' contains non-positive values")

        expansion = np.zeros((1, adata.n_vars))
        for culture, _, c_idx in _anndata_groupby(adata, "__culture_index__"):
            t0_idx = c_idx & t0_mask
            if t0_idx.sum() == 0:
                raise ValueError(f"No t0 sample found for culture '{culture}'")
            growth_t0_mean = growth[t0_idx].mean()
            if growth_type == "density":
                expansion[:, c_idx] = np.log(growth[c_idx]) - np.log(growth_t0_mean)
            else:
                expansion[:, c_idx] = growth[c_idx] * np.log(2.)
        expansion = -expansion
    else:
        raise ValueError(
            "No expansion axis available: set use_spike=True with a spike strain, "
            "or provide growth_column with per-sample density measurements."
        )

    adata.var["__log_expansion__"] = expansion.squeeze(axis=0)
    return adata
