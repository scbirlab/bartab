"""Command-line interface for bartab."""

from argparse import FileType, Namespace
import sys

from carabiner import print_err
from carabiner.cliutils import CLIApp, CLICommand, CLIOption, clicommand

from . import __version__, appname


@clicommand(message="Calculating fitness with the following parameters")
def _fitness(args: Namespace) -> None:

    from .io import load_anndata
    from .models.anndata import AnnDataWLSModel, AnnDataHillModel
    from .transforms import compute_log_ratios

    adata = load_anndata(
        counts=args.inputs,
        sample_meta=args.sample_sheet,
        strain_meta=args.barcode_sheet,
        reference=args.reference,
        count_column=args.count_column,
        timepoint_column=args.timepoint_column,
        t0=args.t0,
        concentration_column=args.concentration_column,
        strain_id=args.barcode_column,
        sample_id=args.sample_column,
        culture_id=args.culture_column,
        spike=args.spike_name,
    )
    print_err(adata)
    adata = compute_log_ratios(
        adata=adata,
        pseudocount=args.pseudocount,
        volume_column=args.volume_column,
        growth_column=args.growth_column,
        growth_type=args.growth_type,
        use_spike=args.use_spike,
    )
    print_err(adata)
    print_err(f"[INFO] Using a weighted least squares model")
    model = AnnDataWLSModel()
    kwargs = {}
    results = model.fit(adata=adata, **kwargs)
    if args.concentration_column is not None:
        if args.concentration_column not in adata.var:
            raise ValueError(
                f"Concentration column '{args.concentration_column}' must be in the --sample-sheet."
            )
        print_err(
            "[INFO] Using a Hill non-linear model with "
            f"'{args.concentration_column}' for concentration."
        )
        model = AnnDataHillModel()
        kwargs = {"concentration": args.concentration_column}
        results = model.fit(adata=results, **kwargs)
    print_err(results)

    print_err(f"[INFO] Writing to {args.output}")
    if args.output.removesuffix(".gz").endswith(".h5ad"):
        results.write(
            args.output,
            compression=(
                "gzip" if args.output.endswith(".gz") 
                else None
            ),
        )
        results.obs.to_csv(
            args.output.removesuffix(".gz").removesuffix(".h5ad") + ".csv",
        )
    elif args.output.removesuffix(".gz").endswith((".csv", ".tsv", ".txt")):
        results.obs.to_csv(
            args.output,
            sep=(
                "," if args.output.removesuffix(".gz").endswith(".csv") 
                else "\t"
            ),
        )
    return None


@clicommand(message="Plotting with the following parameters")
def _plot(args: Namespace) -> None:

    from anndata.io import read_h5ad
    from .plotting import (
        dose_response, 
        expansion_vs_count,
        expansion_vs_ratio, 
        pred_vs_true,
        time_vs_count, 
        time_vs_ratio, 
        volcano
    )

    print_err(
        f"[INFO] Loading from {args.inputs}"
    )
    adata = read_h5ad(args.inputs)
    if args.model_type == "HillFitnessModel":
        fig, axes = dose_response(
            adata,
            highlight_barcodes=args.highlight,
            model_name="WLS",
            filename=args.output + f"_dr.{args.plot_format}",
        )
    fig, axes = time_vs_count(
        adata,
        highlight_barcodes=args.highlight,
        filename=args.output + f"_time-count.{args.plot_format}",
    )
    fig, axes = time_vs_ratio(
        adata,
        highlight_barcodes=args.highlight,
        filename=args.output + f"_time-ratio.{args.plot_format}",
    )
    fig, axes = expansion_vs_count(
        adata,
        highlight_barcodes=args.highlight,
        filename=args.output + f"_expansion-count.{args.plot_format}",
    )
    fig, axes = expansion_vs_ratio(
        adata,
        highlight_barcodes=args.highlight,
        filename=args.output + f"_expansion-ratio.{args.plot_format}",
    )
    fig, axes = pred_vs_true(
        adata,
        highlight_barcodes=args.highlight,
        model_name=args.model_type,
        filename=args.output + f"_pred-obs.{args.plot_format}",
    )

    fig, axes = volcano(
        adata,
        highlight_barcodes=args.highlight,
        model_name=args.model_type,
        param="ic50" if args.model_type == "HillFitnessModel" else "fitness",
        xscale="log" if args.model_type == "HillFitnessModel" else "linear",
        vline=None if args.model_type == "HillFitnessModel" else 1.,
        p="log_ic50_p" if args.model_type == "HillFitnessModel" else "slope_p",
        filename=args.output + f"_volcano.{args.plot_format}",
    )
    return None


@clicommand(message="Simulating with the following parameters")
def _simulate(args: Namespace) -> None:

    import json

    import numpy as np
    import pandas as pd
    from scipy.special import expit as sigmoid

    from .simulation import calculate_growth_curves, reads_sampler

    rng = np.random.default_rng(args.seed)
    print_err(
        f"[INFO] Loading from {args.inputs}"
    )
    with open(args.inputs) as f:
        fitness_input = json.load(f)
    if args.generate_controls:
        fitness_input |= {
            f"ctrl_{i:03d}": [10000., 1.] if args.n_dose else 1.
            for i in range(args.generate_controls)
        }
    if args.generate_more:
        fitness_input |= {
            f"str_{i:03d}": [v, 0.] if args.n_dose else v
            for i, v in enumerate(
                rng.uniform(0., 1.5, size=args.generate_more)
            )
        }
    if args.n_dose:
        doses = np.geomspace(
            args.dose_max / ((args.n_dose - 1) ** args.dose_fold),
            args.dose_max,
            num=args.n_dose,
        )
        fitness = {
            dose: {
                strain: bottom + (1. - bottom) * sigmoid(1. * (np.log(ic50) - np.log(dose)))
                for strain, (ic50, bottom) in fitness_input.items()
            } for dose in doses
        }
    else:
        fitness = {1.: fitness_input}

    count_df = []
    meta_df = []
    strain_meta_df = []
    sample_ids = []
    
    for i, (dose, _fitness) in enumerate(fitness.items()):
        strains = list(_fitness.keys())
        t, ref_expansion, growths = calculate_growth_curves(
            inoculum=args.inoculum, 
            inoculum_var=.1, 
            carrying_capacity=args.carrying_capacity, 
            fitness=_fitness, 
            n_timepoints=args.timepoints,
            max_time=args.max_time,
            seed=args.seed,
        )
        assert len(t) == args.timepoints == len(growths[strains[0]]), print_err(len(t), args.timepoints, len(growths[strains[0]]))
        
        counts = reads_sampler(
            growths, 
            sample_frac=.1, 
            seq_depth=args.reads_per_barcode, 
            reps=args.n_cultures, 
            variance=.001,
            seed=args.seed,
        )
        # print_err(counts.shape)  # strains x reps x time 

        # --- count table: strains × samples ---
        growth_t0 = sum(v[0] for _, v in growths.items())
        for strain_name, strain_arr in zip(_fitness, counts):
            for r, rep_arr in enumerate(strain_arr):
                for j, (_t, time_arr) in enumerate(zip(t, rep_arr)):
                    replicate_id =  f"dose_{i}-rep_{r}"
                    sample_id = f"{replicate_id}-t_{j}"
                    obs_id = f"{sample_id}-strain_{strain_name}"
                    count_df.append({
                        "sample_id": sample_id,
                        "obs_id": obs_id,
                        "strain_id": strain_name,
                        "count": time_arr,
                    })
                    this_sample_growth = sum(v[j] for _, v in growths.items())
                    meta_df.append({
                        "sample_id": sample_id,
                        "replicate": replicate_id,
                        "timepoint": _t,
                        "dose": dose,
                        "growth": this_sample_growth,
                        "generations": np.log(this_sample_growth / growth_t0) / np.log(2.),
                        "is_t0": _t == 0.,
                    })
                    strain_meta_df.append({
                        "strain_id": strain_name,
                    })
    
    counts_df = pd.DataFrame(count_df)
    meta_df = pd.DataFrame(meta_df).drop_duplicates()
    strain_meta_df = pd.DataFrame(strain_meta_df).drop_duplicates()

    strain_meta_df = strain_meta_df.assign(
        **{(
            "dose_params" if args.n_dose else "true_fitness"
            ): lambda x: x["strain_id"].map(fitness_input)
        },
        is_reference=lambda x: x["strain_id"] == "wt",
        is_spike=lambda x: x["strain_id"] == "spike",
    )

    counts_df.to_csv(f"{args.output}_count.csv", index=False)
    meta_df.to_csv(f"{args.output}_sample_meta.csv", index=False)
    strain_meta_df.to_csv(f"{args.output}_strain_meta.csv", index=False)

    print_err("\nCounts:")
    print_err(counts_df)
    print_err("\nSample metadata:")
    print_err(meta_df)
    print_err("\nStrain metadata (ground truth):")
    print_err(strain_meta_df)
    return None


def main() -> None:
    """Main function for CLI app.
    """
    inputs_named = CLIOption(
        'inputs',
        type=str,
        help='Count table(s) in TSV or CSV format.',
    )
    inputs = CLIOption(
        'input', 
        type=FileType('r'),
        default=sys.stdin,
        nargs='?',
        help='Input file. Default: STDIN.',
    )
    outputs_named = CLIOption(
        '--output', '-o', 
        type=str,
        required=True,
        help='Output file.',
    )
    outputs = CLIOption(
        '--output', '-o', 
        type=FileType('w'),
        default=sys.stdout,
        help='Output file. Default: STDOUT',
    )

    formatting = CLIOption(
        '--format', '-f', 
        type=str,
        default='TSV',
        choices=['TSV', 'CSV', 'tsv', 'csv'],
        help='Format of files.',
    )

    conditions = CLIOption(
        '--sample-sheet', '-c', 
        type=str,
        required=True,
        help='Table mapping sample IDs to timepoints and possibly growth and inducer concentration.',
    )
    barcodes = CLIOption(
        '--barcode-sheet', '-b', 
        type=str,
        required=True,
        help='Table mapping strain IDs and any metadata.',
    )
    reference = CLIOption(
        '--reference', '-r',
        type=str,
        required=True,
        help="Non-targeting guide name(s) to use for reference.",
    )
    count_column = CLIOption(
        '--count-column', '-a',
        type=str,
        default="count",
        help='Column of count table corresponding to counts.',
    )
    barcode_column = CLIOption(
        '--barcode-column', '-u',
        type=str,
        default="strain_name",
        help='Column of count table corresponding to barcode/guide/strain names.',
    )
    highlight_barcodes = CLIOption(
        '--highlight', '-X',
        type=str,
        default=None,
        nargs="*",
        help='Names of barcode/guide/strains to highlight in plots.',
    )
    culture_column = CLIOption(
        '--culture-column',
        type=str,
        default="culture_id",
        help='Column of count table corresponding to individual culture repicates.',
    )
    growth_column = CLIOption(
        '--growth-column', '-d',
        type=str,
        default="OD600",
        help='Column of sample table corresponding to growth measurement.',
    )
    growth_type = CLIOption(
        '--growth-type',
        type=str,
        choices=["density", "generations"],
        default="density",
        help='Type of growth in the table under --growth-column.',
    )
    concentration_column = CLIOption(
        '--concentration-column', '-k',
        type=str,
        default=None,
        help='Column of conditions table corresponding to concentration.',
    )
    model_type = CLIOption(
        '--model-type',
        default="WLS",
        choices=["WLS", "OLS", "HillFitnessModel"],
        help='What model type to use.',
    )
    volume_column = CLIOption(
        '--volume-column', 
        type=str,
        default=None,
        help="Column of conditions table corresponding to sample volume. "
        "Used for correcting expansion calculations using spike-in if "
        "volume is not constant between samples.",
    )
    timepoint_column = CLIOption(
        '--timepoint-column', '-t',
        type=str,
        default="timepoint",
        help='Column of conditions table corresponding to timepoint.',
    )
    t0 = CLIOption(
        '--t0', '-0',
        type=float,
        default=None,
        help="Timepoint value of the initial timepoint. "
        "Default: use the minimum recorded timepoint value.",

    )
    sample_column = CLIOption(
        '--sample-column', '-s',
        type=str,
        nargs="*",
        default="sample_id",
        help='Column of conditions table corresponding to sample ID.',
    )
    spike = CLIOption(
        '--spike-name', '-x',
        type=str,
        default=None,
        help="Annotate guide name(s) matching this pattern as spike-ins.",
    )
    use_spike = CLIOption(
        '--use-spike',
        action="store_true",
        help="Instead of growth, use spike-in to estimate culture expansion. "
        "Requires --spike.",
    )
    pseudocount = CLIOption(
        '--pseudocount',
        type=float,
        default=1.,
        help="Value to add to counts before taking a log-transform.",
    )
    controls = CLIOption(
        '--controls', '-z',
        type=str,
        default="controls:",
        help="Control guide name(s) match this pattern.",
    )
    plots = CLIOption(
        '--plot',
        type=str,
        default=None,
        help='If provided, plots will be generated with this filename prefix.',
    )
    plot_format = CLIOption(
        '--plot-format',
        type=str,
        default="png",
        choices=["png", "PNG", "pdf", "PDF"],
        help='File format for plots.',
    )

    ### FOR SIMULATIONS
    generate_controls = CLIOption(
        '--generate-controls', '-z',
        type=int,
        default=0,
        help="Number of control strains (rel. fitness = 1) to generate.",
    )
    generate_more = CLIOption(
        '--generate-more', '-m',
        type=int,
        default=0,
        help="Number of additional strains with random fitness to generate.",
    )
    inoculum = CLIOption(
        '--inoculum', '-n',
        type=int,
        default=1000,
        help="Number of cells per strain in inoculum.",
    )
    carrying_capacity = CLIOption(
        '--carrying-capacity', '-K',
        type=float,
        default=10.,
        help="Maximum fold-change in cell count to support.",
    )
    n_timepoints = CLIOption(
        '--timepoints', '-s',
        type=int,
        default=10,
        help="Number of timepoints to sample.",
    )
    max_time = CLIOption(
        '--max-time', '-t',
        type=int,
        default=10,
        help="Final timepoint to sample.",
    )
    n_cultures = CLIOption(
        '--n-cultures', '-r',
        type=int,
        default=3,
        help="Number of cultures (biological replicates) to generate.",
    )
    seq_depth = CLIOption(
        '--reads-per-barcode', '-b',
        type=int,
        default=1_000,
        help="Mean sequencing depth per barcode per sample.",
    )
    _seed = CLIOption(
        '--seed', '-e',
        type=int,
        default=42,
        help="Seed for reproducible randomness.",
    )
    n_doses = CLIOption(
        '--n-dose', '-d',
        type=int,
        default=None,
        help="Number of doses for inducer dose response, using values in `input` as IC50s. "
        "Default: don't generate dose response.",
    )
    dose_max = CLIOption(
        '--dose-max',
        type=int,
        default=1000,
        help="Maximum dose for inducer dose response.",
    )
    dose_fold = CLIOption(
        '--dose-fold',
        type=int,
        default=2,
        help="Fold dilution for inducer dose response.",
    )
    
    fitness = CLICommand(
        "fit", 
        description="Calculate fitness per strain/barcode.",
        main=_fitness,
        options=[
            inputs_named,
            outputs_named,
            conditions,
            barcodes,
            reference,
            barcode_column,
            culture_column,
            spike,
            count_column,
            sample_column,
            timepoint_column,
            t0,
            growth_column,
            growth_type,
            volume_column,
            use_spike,
            concentration_column,
            model_type,
            pseudocount,
        ],
    )

    plotting = CLICommand(
        "plot", 
        description="Plot modelling results.",
        main=_plot,
        options=[
            inputs_named,
            outputs_named,
            highlight_barcodes,
            model_type,
            plot_format,
        ],
    )

    simulate = CLICommand(
        "sim", 
        description="Simulate read counts from a pooled competition experiment.",
        main=_simulate,
        options=[
            inputs_named,
            outputs_named,
            generate_controls,
            generate_more,
            inoculum,
            carrying_capacity,
            n_timepoints,
            max_time,
            n_doses,
            dose_max,
            dose_fold,
            n_cultures,
            seq_depth,
            _seed,
        ],
    )

    app = CLIApp(
        name=appname,
        version=__version__,
        description="Analysis of barcoded, pooled functional genomics assays.",
        commands=[
            fitness,
            plotting,
            simulate,
        ]
    )
    app.run()
    return None


if __name__ == '__main__':
    main()
