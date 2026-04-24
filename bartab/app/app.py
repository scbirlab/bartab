"""Gradio demo for bartab."""

from typing import Iterable, List, Union
from functools import partial
from io import TextIOWrapper
import os
os.environ["COMMANDLINE_ARGS"] = "--no-gradio-queue"

from carabiner import cast, print_err
from carabiner.decorators import decorator_with_params
from carabiner.pd import read_table
import gradio as gr
import nemony as nm
import numpy as np
import pandas as pd

import anndata
anndata.settings.allow_write_nullable_strings = True
import bartab
from bartab.io import load_anndata
from bartab.models.anndata import AnnDataWLSModel, AnnDataHillModel
from bartab.plotting import (
    dose_response,
    expansion_vs_count,
    expansion_vs_ratio, 
    pred_vs_true,
    time_vs_count, 
    time_vs_ratio, 
    volcano
)
from bartab.transforms import compute_log_ratios

pd.options.future.infer_string = True

MODES: dict = {
    "single": "👟 Fitness in a single condition", 
    "dose response": "📉 Fitness dose response to CRISPRi inducer",
}
def _message(s: str):
    print_err(f"[INFO] {s}")
    gr.Info(s, duration=10)
    return None


def load_input_data(
    filename: str, 
    cols: Iterable
) -> List[pd.DataFrame]:
    df = read_table(filename)
    print_err(df)
    out = [gr.update(value=df, visible=False)]
    for key, col in cols.items():
        if isinstance(col, tuple):
            col_type = col[1]
            if col_type == "string":
                choices = list(df.select_dtypes(include="str"))
            elif col_type == "numeric":
                choices = list(df.select_dtypes(include="number"))
            else:
                choices = list(df)
        else:
            choices = list(df)
        choices = [""] + choices
        print_err(key, f"{choices=}")
        out.append(
            gr.update(
                choices=choices,
                value=key if key in choices else choices[0],
                interactive=True, 
                visible=True,
            )
        )
    print_err(out)
    return out


def load_barcode_names(
    df: pd.DataFrame,
    strain_col: str
) -> List[List[str]]:
    strains = sorted(df[strain_col].unique())
    print_err(strain_col, f"{strains=}")
    return gr.update(
        choices=strains,
        value="wt" if "wt" in strains else strains[0],
        interactive=True, 
        visible=True,
    ), gr.update(
        choices=strains, 
        value="spike" if "spike" in strains else "",
        interactive=True, 
        visible=True,
    ), gr.update(
        choices=strains, 
        value=[],
        allow_custom_value=True,
        interactive=True, 
        visible=True,
    )
    

def _prepare_to_fit(
    counts: pd.DataFrame,
    strain_sheet: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    count_column: str,
    strain_id_column: str,
    timepoint_column: str,
    concentration_column: str,
    sample_id_column: str,
    culture_id_column: str,
    volume_column: str = "volume",
    growth_column: str = "growth",
    reference: str = "wt",
    spike_name: str = "spike",
    spike_mode: str = "Spike",
    growth_type: str = "density",
    pseudocount: float = 1.
):
    use_spike = (spike_mode == "Spike")
    adata = load_anndata(
        counts=counts,
        sample_meta=sample_sheet,
        strain_meta=strain_sheet,
        reference=reference,
        count_column=count_column,
        timepoint_column=timepoint_column,
        # t0=args.t0,
        concentration_column=concentration_column,
        strain_id=strain_id_column,
        sample_id=sample_id_column,
        culture_id=culture_id_column,
        spike=spike_name if spike_name else None,
    )
    print_err(adata)
    adata = compute_log_ratios(
        adata=adata,
        pseudocount=pseudocount,
        volume_column=volume_column,
        growth_column=growth_column,
        growth_type=growth_type,
        use_spike=use_spike,
    )
    print_err(adata)
    return adata


def do_analysis(*args):
    args = [
        a if not (isinstance(a, str) and a == "") else None 
        for a in args
    ]
    mode = args[-1]
    concentration_column = args[6]
    print_err(args[:-1])
    try:
        adata = _prepare_to_fit(*args[:-1])
    except TypeError as e:
        print_err(*args[:-1])
        raise e
    _message("Using a weighted least squares model")
    model = AnnDataWLSModel()
    results = model.fit(adata=adata)
    if mode == MODES["dose response"]:
        _message(
            "Using a Hill non-linear model with "
            f"'{concentration_column}' for concentration."
        )
        model = AnnDataHillModel()
        results = model.fit(adata=results, concentration=concentration_column)
    return gr.update(value=results.obs, visible=True), results


def _fig2img(fig):
    import PIL
    # img = PIL.Image.frombytes(
    #     "RGBa", 
    #     fig.canvas.get_width_height(),
    #     fig.canvas.buffer_rgba(),
    # )
    import io
    buf = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = PIL.Image.open(buf)
    return img


@decorator_with_params
def _plot_wrapper(fn, message="Plotting..."):
    def _fn(*args, **kwargs):
        if args[1] == "":
            args[1] = None
        if message:
            _message(message)
        fig, axes = fn(*args, **kwargs)
        if isinstance(fig, tuple) and isinstance(axes, bool):
            fig, vis = fig
            return gr.update(
                value=_fig2img(fig) if fig is not None else fig, 
                visible=vis,
            )
        else:
            return gr.update(value=_fig2img(fig), visible=True)
    return _fn


@_plot_wrapper()
def _plot_dose_response(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    mode: str = MODES["single"]
):
    do_dose_response = mode == MODES["dose response"]
    if do_dose_response:
        print_err("Plotting dose response")
        fig, axes = dose_response(
            adata,
            highlight_barcodes=highlight,
            model_name="WLS",
            control_prefix=control_prefix,
        )
        return fig, axes
    else:
        print_err("Skipping dose response")
        return (None, None), False


@_plot_wrapper(message="Plotting time vs count")
def _plot_time_vs_count(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    *args, **kwargs
):
    return time_vs_count(
        adata,
        highlight_barcodes=highlight,
        control_prefix=control_prefix,
    )


@_plot_wrapper(message="Plotting expansion vs count")
def _plot_expansion_vs_count(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    *args, **kwargs
):
    return expansion_vs_count(
        adata,
        highlight_barcodes=highlight,
        control_prefix=control_prefix,
    )


@_plot_wrapper(message="Plotting time vs ratio")
def _plot_time_vs_ratio(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    *args, **kwargs
):
    return time_vs_ratio(
        adata,
        highlight_barcodes=highlight,
        control_prefix=control_prefix,
    )


@_plot_wrapper(message="Plotting expansion vs ratio")
def _plot_expansion_vs_ratio(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    *args, **kwargs
):
    return expansion_vs_ratio(
        adata,
        highlight_barcodes=highlight,
        control_prefix=control_prefix,
    )

@_plot_wrapper(message="Plotting predicted vs observed")
def _plot_pred_vs_true(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    mode: str = MODES["single"]
):
    do_dose_response = mode == MODES["dose response"]
    return pred_vs_true(
        adata,
        model_name="HillFitnessModel" if do_dose_response else "WLS",
        highlight_barcodes=highlight,
        control_prefix=control_prefix,
    )


@_plot_wrapper(message="Plotting volcano")
def _plot_volcano(
    adata,
    highlight=None,
    control_prefix: str = "ctrl_",
    mode: str = MODES["single"]
):
    do_dose_response = mode == MODES["dose response"]
    return volcano(
        adata,
        highlight_barcodes=highlight,
        control_prefix=control_prefix,
        model_name="HillFitnessModel" if do_dose_response else "WLS",
        param="ic50" if do_dose_response else "fitness",
        xscale="log" if do_dose_response else "linear",
        vline=None if do_dose_response else 1.,
        p="log_ic50_p" if do_dose_response else "slope_p",
    )


def download_tables(
    df: pd.DataFrame,
    adata
) -> str:
    df_hash = nm.hash(pd.util.hash_pandas_object(df).values)
    filename = f"bartab-{df_hash}"
    filename_csv = f"{filename}.csv"
    df.to_csv(filename, index=False)
    filename_adata = f"{filename}.h5ad"
    adata.write(filename_adata)
    return gr.update(
        value=filename_csv, 
        visible=True,
    ), gr.update(
        value=filename_adata, 
        visible=True,
    )



def _file_input(**kwargs):
    return partial(gr.File,
        file_types=[".xlsx", ".csv", ".tsv", ".txt"],
    )(**kwargs)


def _load_from_file(*args):
    if len(args) > 1:
        return args
    else:
        return args[0]


def _invisible_dropdown(**kwargs):
    return partial(gr.Dropdown,
        choices=[],
        interactive=False,
        visible=True,
    )(**kwargs)


def _invisible_plot(**kwargs):
    return partial(gr.Image,
        visible=False,
    )(**kwargs)


with gr.Blocks() as demo:
    gr.Markdown(
        f"""
        # 🍹 bartab: Fitness from pooled competition assays

        *Using* `bartab` v{bartab.__version__} | [Documentation](https://github.com/scbirlab/bartab) | [Tutorial on analysis principles](https://huggingface.co/spaces/scbirlab/tutorial-seq-fitness)
        
        Infer the competitive fitness of barcoded strains from next-generation 
        sequencing of pooled growth experiments.
        
        Upload your count table, sample sheet, and barcode sheet, then 
        click **Calculate fitness**.

        """
    )
    gr.Markdown(
        """
        ---

        ## 1️⃣ Input tables

        Three tables are required. You can upload CSV, TSV, or XLSX files, 
        or try one of the **example datasets** below.

        """
    )
    input_filenames = {
        "count_table": gr.Textbox(interactive=False, visible=False),
        "sample_sheet": gr.Textbox(interactive=False, visible=False),
        "strain_sheet": gr.Textbox(interactive=False, visible=False),
    }
    app_root = os.path.dirname(__file__)
    data_path = os.path.join(app_root, "data", "examples", "single-point")
    control_strains = {
        "reference": _invisible_dropdown(
            label="Reference (WT) barcode name",
            render=False,
        ),
        "spike": _invisible_dropdown(
            label="Spike-in barcode name (if using)",
            render=False,
        ),
    }
    analysis_opts = {
        "use_spike": gr.Radio(
            label="Culture expansion uses:",
            choices=["Spike", "Growth"],
            value="Spike",
            render=False,
        ),
        "growth_type": gr.Radio(
            label="Growth type",
            choices=["density", "generations"],
            value="density",
            visible=False,
            render=False,
        ),
        "pseudocount": gr.Number(
            label="Pseudocount",
            value=1.,
            render=False,
        ),
    }
    plotting_opts = {
        "highlight": _invisible_dropdown(
            label="Strain(s) to highlight in plots",
            render=False,
            multiselect=True,
        ),
        "controls": gr.Textbox(
            label="Prefix of control barcode names",
            value="ctrl_",
            render=False,
        ),
    }
    mode_switch = gr.Radio(
        label="Analysis mode",
        choices=list(MODES.values()),
        value=MODES["single"],
        render=False,
    )
    examples = gr.Examples(
        label="Examples with synthetic data",
        examples=[
            [
                os.path.join(data_path, "test_count.csv"),
                os.path.join(data_path, "test_sample_meta.csv"),
                os.path.join(data_path, "test_strain_meta.csv"),
                MODES["single"],
                "Spike",
            ],
            [
                os.path.join(data_path, "test_count.csv"),
                os.path.join(data_path, "test_sample_meta.csv"),
                os.path.join(data_path, "test_strain_meta.csv"),
                MODES["single"],
                "Growth",
            ],
            [
                os.path.join(data_path, "dose-response_count.csv"),
                os.path.join(data_path, "dose-response_sample_meta.csv"),
                os.path.join(data_path, "dose-response_strain_meta.csv"),
                MODES["dose response"],
                "Spike",
            ],
            [
                os.path.join(data_path, "dose-response_count.csv"),
                os.path.join(data_path, "dose-response_sample_meta.csv"),
                os.path.join(data_path, "dose-response_strain_meta.csv"),
                MODES["dose response"],
                "Growth",
            ],
        ],
        example_labels=[
            ["Single point, using spike-in"],
            ["Single point, using growth"],
            ["Dose response, using spike-in"],
            ["Dose response, using growth"],
        ],
        inputs=[
            input_filenames["count_table"], 
            input_filenames["sample_sheet"], 
            input_filenames["strain_sheet"],
            mode_switch,
            analysis_opts["use_spike"],
        ],
        # cache_examples=True,
        # cache_mode="eager",
    )
    input_files = {}
    input_cols = {}
    go_button = gr.Button(
        value="🚀 Calculate fitness!",
        interactive=False,
        render=False,
    )
    with gr.Row():
        with gr.Column():
            gr.Markdown(
                """
                ---

                ### 🧮 Count table

                One row per barcode per sample. Must contain:
                - a column of **barcode/strain identifiers** (matching your barcode sheet)
                - a column of **sample identifiers** (matching your sample sheet)  
                - a column of **read or UMI counts**

                """
            )
            input_files["count_table"] = _file_input(
                label="Upload your barcode sequencing counts data here",
            )
            input_cols["count_table"] = {
                "count": (_invisible_dropdown(
                    label="Counts column",
                ), "numeric"),
            }
        with gr.Column():
            gr.Markdown(
                """
                ---

                ### 📶 Barcode information

                One row per unique barcode. Must contain:
                - a column of **barcode identifiers**

                Optionally: any metadata about strains (gene targets, constructs, 
                etc.). These will be carried through to the output.

                """
            )
            input_files["strain_sheet"] = _file_input(
                label="Upload your barcode information here",
            )
            input_cols["strain_sheet"] = {
                "strain_id": (_invisible_dropdown(
                    label="Barcode identifier column",
                ), "string"),
            }
    with gr.Row():
        gr.Markdown(
            r"""
            ---

            ### 🧪 Sample sheet

            One row per sample (sequencing library). Must contain:
            - **Sample ID**: unique identifier matching the count table
            - **Culture ID**: biological replicate identifier. Samples from 
            the same culture share this label
            - **Timepoint**: numeric timepoint values. The earliest timepoint 
            is treated as $t_0$.

            **For spike-in normalisation**: 
            no extra columns needed. Just include your spike-in barcode in 
            the count table and barcode sheet.

            **For growth-based normalisation**: add a column of OD600, CFU/mL, 
            or generation counts measured at each sample.

            **For dose-response analysis**: add a column of inducer/drug 
            concentrations. Samples with concentration = 0 are treated as 
            uninduced controls.

            **For adaptive-volume sampling** (if you took different volumes 
            from each sample): add a column of sampled volumes.

            """
        )
        input_files["sample_sheet"] = _file_input(
            label="Upload your sample information here",
        )
    with gr.Row():
        input_cols["sample_sheet"] = {
            "timepoint": (_invisible_dropdown(
                label="Timepoint column",
            ), "any"),
            "dose": (_invisible_dropdown(
                label="Concentration column",
            ), "numeric"),
            "sample_id": (_invisible_dropdown(
                label="Individual sample ID column",
            ), "string"),
            "replicate": (_invisible_dropdown(
                label="Culture / biological replicate column",
            ), "string"),
            "volume": (_invisible_dropdown(
                label="Volume column (if using)",
            ), "numeric"),
            "growth": (_invisible_dropdown(
                label="Growth column (if using)",
            ), "numeric"),            
        }
        
    with gr.Row():
        input_data = {
            key: gr.Dataframe(
                label=f"Input data: {key}",
                max_height=50,
                visible=False,
                interactive=False,
            ) for key in input_files
        }

    adata = gr.State()
    with gr.Row():
        with gr.Column():
            gr.Markdown(
                """
                ---

                ## 2️⃣ Control strains

                - **Reference (WT)**: the strain relative to which all 
                fitness values are calculated. Fitness = 1 by definition.
                - **Spike-in**: a non-growing strain (e.g. heat-killed or 
                plasmid-only) added at a fixed concentration before library 
                preparation. Used to infer how much the reference strain has 
                expanded between timepoints, removing the need for growth 
                measurements. Leave blank if using growth measurements instead.

                """
            )
            for key, val in control_strains.items():
                val.render()
        with gr.Column():
            gr.Markdown(
                """
                ---

                ## 3️⃣ Analysis options

                - **Culture expansion**: choose **Spike** if you have a 
                non-growing spike-in control, or **Growth** if you have 
                OD600/CFU measurements.
                - **Pseudocount**: added to all counts before log transformation 
                to avoid log(0).
                - **Analysis mode**: choose **Single concentration** for standard fitness 
                screens, or **Dose response** if your sample sheet contains a concentration 
                column. Dose response fitting uses a 2-parameter Hill model to estimate the 
                IC₅₀ and maximum effectfor each barcode.

                """
            )
            for key, val in analysis_opts.items():
                val.render()
        with gr.Column():
            gr.Markdown(
                """
                ---
                
                ## Plotting options

                """
            )
            for key, val in plotting_opts.items():
                val.render()

    with gr.Column():
        mode_switch.render()
        go_button.render()
        
    mode_switch.change(
        lambda x: gr.update(value=x),
        inputs=[mode_switch],
        outputs=[go_button],
    )
    gr.Markdown(
        r"""
        ## 4️⃣ Results

        Fitness values are estimated by weighted least squares regression of the log-ratio of each barcode against the reference strain, using the spike-in or growth measurements as the x-axis. 

        **Key output columns**:
        - For single concentration:
            - `fitness`: relative fitness ($w_i / w_{wt}$). Values < 1 indicate growth disadvantage; > 1 indicates advantage.
            - `fitness_low` / `fitness_high`: 95% confidence interval bounds.
            - `slope_p`: p-value for the slope being different from 0 (i.e. fitness ≠ 1).
        - For dose-response: 
            - `log_ic50` (log₁₀ concentration at 50% inhibition) and `log_ic50_p`.

        Results and the full annotated dataset (`.h5ad`) can be downloaded below.

        """
    )

    plots = {}
    with gr.Row():
        plots |= {
            "dose_response": (
                _invisible_plot(label="Dose response"),
                _plot_dose_response,
            ),
            "count_time": (
                _invisible_plot(label="Time vs count"),
                _plot_time_vs_count,
            ),
            "count_exp": (
                _invisible_plot(label="Expansion vs count"),
                _plot_expansion_vs_count,
            ),
        }
    with gr.Row():
        plots |= {
            "count_exp": (
                _invisible_plot(label="Time vs ratio"),
                _plot_time_vs_ratio,
            ),
            "ratio_exp": (
                _invisible_plot(label="Expansion vs ratio"),
                _plot_expansion_vs_ratio,
            ),
        }
    with gr.Row():
        plots |= {
            "pred_obs": (
                _invisible_plot(label="Predicted vs observed"),
                _plot_pred_vs_true,
            ),
            "volcano": (
                _invisible_plot(label="Volcano"),
                _plot_volcano,
            ),
        }
    with gr.Row():
        download = gr.DownloadButton(
            label="Download parameters as CSV",
            visible=False,
        )
        download_adata = gr.DownloadButton(
            label="Download all analysis as .h5ad",
            visible=False,
        )
    output_table = gr.Dataframe(
        label="Fitted parameters",
        # max_height=100,
        visible=False,
        interactive=False,
    )

    # ======
    # EVENTS
    # ======

    for key, input_file in input_files.items():
        input_columns = input_cols[key]
        event_fn = {
            "fn": partial(load_input_data, cols=input_columns),
            "outputs": [input_data[key]] + [
                col[0] if isinstance(col, tuple) else col 
                for _, col in input_columns.items()
            ],
        }
        input_filenames[key].change(
            _load_from_file, 
            inputs=[input_filenames[key]], 
            outputs=[input_file],
        ).then(
            **event_fn, 
            inputs=[input_filenames[key]], 
        )
        input_file.upload(
            **event_fn, 
            inputs=[input_file], 
        )
        
    input_cols["strain_sheet"]["strain_id"][0].change(
        load_barcode_names,
        inputs=[
            input_data["strain_sheet"], 
            input_cols["strain_sheet"]["strain_id"][0],
        ],
        outputs=[
            control_strains["reference"],
            control_strains["spike"],
            plotting_opts["highlight"],
        ],
    ).then(
        lambda : gr.update(interactive=True),
        inputs=[],
        outputs=[go_button],
    )

    analysis_opts["use_spike"].change(
        lambda x: gr.update(visible=x == "Growth"),
        inputs=[analysis_opts["use_spike"]],
        outputs=[analysis_opts["growth_type"]],
    )

    comptuation_inputs = (
        [
            f for _, f in input_data.items()
        ] + [
            opt[0] if isinstance(opt, tuple) else opt
            for _, input_col in input_cols.items()
            for _, opt in input_col.items()
        ] + [
            v for k, v in control_strains.items()
        ] + [
            v for k, v in analysis_opts.items()
        ]
    )

    evt = go_button.click(
        fn=do_analysis,
        inputs=comptuation_inputs + [mode_switch],
        outputs=[
            output_table,
            adata,
        ],
    ).then(
        download_tables,
        inputs=[output_table, adata],
        outputs=[download, download_adata],
    )
    
    for key, (p, fn) in  plots.items():
        evt = evt.then(
            fn,
            inputs=[
                adata, 
                plotting_opts["highlight"], 
                plotting_opts["controls"], 
                mode_switch,
            ],
            outputs=[p],
        )
    
if __name__ == "__main__":
    demo.queue() 
    demo.launch(share=True)
