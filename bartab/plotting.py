"""Plotting functions."""

from typing import Callable
import os

from carabiner.mpl import set_plot_palette, grid, figsaver, scattergrid
import numpy as np
import pandas as pd


def save_plot(
    fig,
    filename,
    **kwargs
):
    path, ext = os.path.splitext(filename)
    return figsaver(
        format=ext.lstrip("."), 
        dpi=600,
    )(fig, path, **kwargs)


def scatter(
    ax,
    x: str,
    y: str,
    data,
    scatter_opts: dict = None,
    vline: float = None,
    hline: float = None,
    xlabel: str =None,
    ylabel: str =None,
    **kwargs
):
    defaults = {
        "facecolor": "none",
        "edgecolor": "C0",
        "linewidth": .5,
        "s": 10.,
        "label": "_none",
    }
    if scatter_opts is None:
        scatter_opts = {}
    ax.scatter(
        data[x],
        data[y],
        **(defaults | scatter_opts),
    )
    if vline:
        ax.axvline(vline, color="lightgrey", zorder=-1)
    if hline:
        ax.axhline(hline, color="lightgrey", zorder=-1)
    
    ax.set(
        xlabel=xlabel or x,
        ylabel=ylabel or x,
        **kwargs,
    )
    return ax


def volcano(
    adata,
    model_name,
    param: str = "fitness",
    p: str = "slope_p",
    filename: str = None
):
    fig, ax = grid()
    df = adata.obs
    x = f"{model_name}:{param}"
    y = f"{model_name}:{p}"
    ax = scatter(
        ax=ax,
        x=x,
        y=y,
        data=df,
    )
    ax = scatter(
        ax=ax,
        x=x,
        y=y,
        data=df.query("strain_id.str.startswith('ctrl_')"),
        scatter_opts={
            "facecolor": "lightgrey",
            "edgecolor": "none",
        },
        vline=1.,
        hline=.05,
        xlabel="Relative fitness",
        ylabel="P",
        # xscale="log",
        # yscale="log",
    )
    ax = scatter(
        ax=ax,
        x=x,
        y=y,
        data=df.query("__is_spike__"),
        scatter_opts={
            "facecolor": "C1",
            "edgecolor": "none",
        },
        vline=1.,
        hline=.05,
        xlabel=f"Parameter: {param}",
        ylabel=f"P: {p}",
        # xscale="log",
        # yscale="log",
    )
    if filename is not None:
        save_plot(fig, filename, df=df[[x, y]].reset_index())
    return fig, ax


def _layered_scatter_barcodes(
    adata, 
    layer: str, 
    x_obs: str = None,
    x_layer: str = None,
    filename: str = None,
    palette: str = None,
    spike_label: str = "__is_spike__",
    control_prefix="ctrl_",
    color_by_barcode: bool = False,
    default_color: str = "lightgrey",
    exp_x: bool = False, 
    exp_y: bool = False,
    callback=None,
    **kwargs
):
    set_plot_palette("vibrant")
    if palette is not None:
        set_plot_palette(palette)
    fig, ax = grid()

    dfs = []
    if x_obs is not None:
        _x = adata.var[x_obs].values
    elif x_layer is None:
        _x = adata.X
    else:
        _x = adata.layers[x_layer]
    if layer is None:
        Y = adata.X
    else:
        Y = adata.layers[layer]
    # print(_x.shape, Y.shape)
    X = np.broadcast_to(_x[None] if _x.ndim < Y.ndim else _x, Y.shape).copy()
    for i, (_x, _y, idx) in enumerate(zip(
        X,
        Y, 
        adata.obs.index
    )):
        if idx.startswith(control_prefix):
            scatter_opts = {
                "facecolor": "dimgrey",
                "edgecolor": "none",
                "zorder": 3,
                "s": 15.,
            }
        elif adata.obs.loc[idx][spike_label]:
            scatter_opts = {
                "facecolor": "C1",
                "edgecolor": "none",
                "zorder": 10,
                "s": 15.,
            }
        else:
            scatter_opts = {
                "edgecolor": f"C{i}" if color_by_barcode else default_color,
                "s": 5.,
            }
        this_df = pd.DataFrame({
            x_obs: np.exp(_x) if exp_x else _x, 
            layer: np.exp(_y) if exp_y else _y,
        }).assign(**{
            "__index__": idx, 
            "__i__": i, 
        })
        ax = scatter(
            ax=ax,
            x=x_obs,
            y=layer,
            data=this_df,
            scatter_opts=scatter_opts,
        )
        dfs.append(this_df)
    dfs = pd.concat(dfs)
    if callback is not None and isinstance(callback, Callable):
        callback(ax)
    ax.set(**({
        "xscale": "log" if exp_x else "linear", 
        "yscale": "log" if exp_y else "linear",
    } | kwargs))
    if filename is not None:
        save_plot(fig, filename, df=dfs)
    return fig, ax


def pred_vs_true(
    adata,
    model_name,
    filename: str = None,
    **kwargs
):
    fig, ax = _layered_scatter_barcodes(
        adata, 
        x_layer=f"{model_name}:predicted", 
        layer=None, 
        filename=filename,
        exp_x=False, 
        exp_y=False,
        callback=lambda ax: ax.plot(ax.get_xlim(), ax.get_xlim(), color="lightgrey", zorder=-1),
        xlabel=f"Predicted: {model_name}",
        ylabel="Observed",
        xscale="log",
        yscale="log",
        **kwargs,
    )
    return fig, ax


def expansion_vs_ratio(
    adata,
    filename: str = None,
    **kwargs
):
    fig, ax = _layered_scatter_barcodes(
        adata, 
        layer="__log_ratio__", 
        x_obs="__log_expansion__", 
        filename=filename,
        exp_x=True, 
        exp_y=True,
        callback=lambda ax: ax.axhline(0., color="lightgrey", zorder=-1),
        xlabel="Culture expansion",
        ylabel="Barcode expansion vs ref",
        **kwargs,
    )
    return fig, ax


def time_vs_ratio(
    adata,
    filename: str = None,
    **kwargs
):
    fig, ax = _layered_scatter_barcodes(
        adata, 
        layer="__log_ratio__", 
        x_obs="timepoint", 
        filename=filename,
        exp_x=False, 
        exp_y=True,
        callback=lambda ax: ax.axhline(0., color="lightgrey", zorder=-1),
        xlabel="Timepoint",
        ylabel="Barcode expansion vs ref",
        **kwargs,
    )
    return fig, ax


def time_vs_count(
    adata,
    filename: str = None,
    **kwargs
):
    fig, ax = _layered_scatter_barcodes(
        adata, 
        layer=None, 
        x_obs="timepoint", 
        filename=filename,
        exp_x=False, 
        exp_y=False,
        callback=lambda ax: ax.axhline(0., color="lightgrey", zorder=-1),
        xlabel="Timepoint",
        ylabel="Count",
        yscale="log",
        **kwargs,
    )
    return fig, ax

