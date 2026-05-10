"""Plotting functions."""

from typing import Callable, Iterable, Optional
import os

from carabiner.mpl import add_legend, set_plot_palette, grid, figsaver, scattergrid
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
        "edgecolor": "lightgrey",
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
    if "label" in scatter_opts:
        if all([
            len(scatter_opts["label"]) > 0,
            scatter_opts["label"][0] != "_",
        ]):
            add_legend(ax)
    return ax


def _avoid_color_collision(i, avoid=None):
    if avoid is None:
        avoid = []
    if i in avoid:
        return _avoid_color_collision(i + 1, avoid)
    else:
        return i


def volcano(
    adata,
    model_name,
    param: str = "fitness",
    p: str = "slope_p",
    control_prefix: str = "ctrl_",
    spike_color: str = "C1",
    highlight_barcodes: Optional[Iterable[str]] = None,
    filename: str = None,
    **kwargs,
):
    fig, ax = grid(aspect_ratio=1.35)
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
        data=df.query("__row_index__.str.startswith(@control_prefix)"),
        scatter_opts={
            "facecolor": "dimgrey",
            "edgecolor": "none",
            "s": 15.,
            "label": "Control",
        },
    )
    ax = scatter(
        ax=ax,
        x=x,
        y=y,
        data=df.query("__is_spike__"),
        scatter_opts={
            "facecolor": spike_color,
            "edgecolor": "none",
            "s": 15.,
            "label": "Spike",
        },
        hline=.05,
        xlabel=f"Parameter: {param}",
        ylabel=f"P: {p}",
        # xscale="log",
        yscale="log",
        **kwargs,
    )
    if highlight_barcodes is None:
        highlight_barcodes = []
    for i, (idx, bc_df) in enumerate(df.query("__row_index__.isin(@highlight_barcodes)").groupby("__row_index__")):
        color = f"C{_avoid_color_collision(i + 2, [int(spike_color[1:]), 7, 8])}"
        ax = scatter(
            ax=ax,
            x=x,
            y=y,
            data=bc_df,
            scatter_opts={
                "facecolor": "none",
                "edgecolor": color,
                "s": 5.,
                "linewidth": 1.,
                "label": idx,
            },
            hline=.05,
            xlabel=f"Parameter: {param}",
            ylabel=f"P: {p}",
            # xscale="log",
            yscale="log",
            **kwargs,
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
    control_prefix: str = "ctrl_",
    color_by_barcode: bool = False,
    default_color: str = "lightgrey",
    spike_color: str = "C1",
    highlight_barcodes: Optional[Iterable[str]] = None,
    exp_x: bool = False, 
    exp_y: bool = False,
    callback=None,
    **kwargs
):
    if highlight_barcodes is None:
        highlight_barcodes = []
    set_plot_palette("vibrant")
    if palette is not None:
        set_plot_palette(palette)
    fig, ax = grid(aspect_ratio=1.35)

    dfs = []
    if x_obs is not None:
        X = adata.var[x_obs].values
        x_col = x_obs
    elif x_layer is None:
        X = adata.X
        x_col = "X" 
    else:
        X = adata.layers[x_layer]
        x_col = x_layer
    if layer is None:
        Y = adata.X
        y_col = "Y" 
    else:
        Y = adata.layers[layer]
        y_col = layer
    # print(_x.shape, Y.shape)
    X = X[None] if X.ndim < Y.ndim else X
    X = np.broadcast_to(X, Y.shape).copy()
    control_plotted, spike_plotted = False, False
    for i, (_x, _y, idx) in enumerate(zip(
        X,
        Y, 
        adata.obs.index,
    )):
        if idx.startswith(control_prefix):
            scatter_opts = {
                "facecolor": "dimgrey",
                "edgecolor": "none",
                "zorder": 3,
                "s": 15.,
                "label": "Controls" if not control_plotted else "_none",
            }
            control_plotted = True
        elif adata.obs.loc[idx][spike_label]:
            scatter_opts = {
                "facecolor": spike_color,
                "edgecolor": "none",
                "zorder": 10,
                "s": 15.,
                "label": "Spike" if not spike_plotted else "_none",
            }
            spike_plotted = True
        elif adata.obs.loc[idx]["__is_reference__"]:
            scatter_opts = {
                "facecolor": "none",
                "edgecolor": "black",
                "zorder": 12,
                "s": 20.,
                "label": "Ref.",
            }
        elif idx in highlight_barcodes:
            color = f"C{_avoid_color_collision(i + 2, [int(spike_color[1:]), 7, 8])}"
            scatter_opts = {
                "edgecolor": color,
                "facecolor": "none",
                "zorder": 10,
                "s": 5.,
                "linewidth": 1.,
                "label": idx,
            }
        else:
            scatter_opts = {
                "edgecolor": f"C{i}" if color_by_barcode else default_color,
                "s": 5.,
                "label": "_none",
            }
        this_df = pd.DataFrame({
            x_col: np.exp(_x) if exp_x else _x, 
            y_col: np.exp(_y) if exp_y else _y,
        }).assign(**{
            "__index__": idx, 
            "__i__": i, 
        })
        ax = scatter(
            ax=ax,
            x=x_col,
            y=y_col,
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
        layer="__log_ratio__", 
        filename=filename,
        exp_x=True, 
        exp_y=True,
        callback=lambda ax: ax.plot(ax.get_xlim(), ax.get_xlim(), color="lightgrey", zorder=-1),
        xlabel=f"Predicted: {model_name}",
        ylabel="Observed:\nbarcode expansion ratio",
        # xscale="log",
        **kwargs,
    )
    return fig, ax



def expansion_vs_count(
    adata,
    filename: str = None,
    **kwargs
):
    fig, ax = _layered_scatter_barcodes(
        adata, 
        layer=None, 
        x_obs="__log_expansion__", 
        filename=filename,
        exp_x=True, 
        exp_y=False,
        # callback=lambda ax: ax.axhline(1., color="lightgrey", zorder=-1),
        xlabel="Culture expansion ratio",
        ylabel="Count",
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
        callback=lambda ax: ax.axhline(1., color="lightgrey", zorder=-1),
        xlabel="Culture expansion ratio",
        ylabel="Barcode expansion ratio",
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
        x_obs=adata.uns["timepoint_column"], 
        filename=filename,
        exp_x=False, 
        exp_y=True,
        callback=lambda ax: ax.axhline(1., color="lightgrey", zorder=-1),
        xlabel="Timepoint",
        ylabel="Barcode expansion ratio",
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
        x_obs=adata.uns["timepoint_column"], 
        filename=filename,
        exp_x=False, 
        exp_y=False,
        # callback=lambda ax: ax.axhline(0., color="lightgrey", zorder=-1),
        xlabel="Timepoint",
        ylabel="Count",
        yscale="log",
        **kwargs,
    )
    return fig, ax


def dose_response(
    adata,
    model_name: str,
    filename: str = None,
    control_prefix: str = "ctrl_",
    highlight_barcodes: Optional[Iterable[str]] = None,
    **kwargs
):
    if highlight_barcodes is None:
        highlight_barcodes = []
    prefix = f"{model_name}:fitness"
    fitnesses = [col for col in adata.obs.columns if col.startswith(prefix)]
    concs = set([float(col.split("@")[-1]) for col in fitnesses])
    fitness_df = (
        adata.obs
        [fitnesses]
        .reset_index()
        .melt(
            id_vars="__row_index__",
            var_name="_param_name",
            value_name="param_value",
        )
        .assign(
            param_name=lambda x: x["_param_name"].str.split(f"{model_name}:").str[1].str.split("@").str[0],
            _concentration=lambda x: x["_param_name"].str.split("@").str[-1].astype(float),
        )
        .pivot(
            index=["__row_index__", "_concentration"],
            columns="param_name",
            values="param_value",
        )
        .reset_index()
    )
    fig, ax = grid(aspect_ratio=1.35)
    control_plotted, spike_plotted = False, False
    for i, (barcode_name, bc_df) in enumerate(
        fitness_df
        .sort_values("_concentration")
        .groupby("__row_index__")
    ):
        if barcode_name == adata.uns["spike"]:
            line_args = {
                "color": "C1",
                "linewidth": 1.,
                "label": "Spike" if not spike_plotted else "_none",
                "zorder": 10,
            }
            spike_plotted = True
        elif barcode_name.startswith(control_prefix):
            line_args = {
                "color": "dimgrey",
                "linewidth": 1.,
                "label": "Controls" if not control_plotted else "_none",
                "zorder": 3,
            }
            control_plotted = True
        elif barcode_name == adata.uns["reference"]:
            scatter_opts = {
                "facecolor": "none",
                "edgecolor": "black",
                "zorder": 12,
                "s": 20.,
                "label": "Ref.",
            }
        elif barcode_name in highlight_barcodes:
            color = f"C{_avoid_color_collision(i + 2, [1, 7, 8])}"
            line_args = {
                "color": color,
                "linewidth": 1.,
                "label": barcode_name,
                "zorder": 10,
            }
        else:
            line_args = {
                "color": "lightgrey",
                "alpha": .4,
                "linewidth": .5,
                "label": "_none",
                "zorder": 0,
            }
        # print(bc_df)
        ax.plot(
            "_concentration",
            "fitness",
            data=bc_df,
            **line_args,
        )
    ax.set(
        xlabel="Concentration",
        ylabel="Fitness",
        xscale="log",
    )
    add_legend(ax)
    if filename is not None:
        save_plot(fig, filename, df=fitness_df)
    return fig, ax

