"""Command-line interface for bartab."""

from argparse import FileType, Namespace
import sys

from carabiner import print_err
from carabiner.cliutils import CLIApp, CLICommand, CLIOption, clicommand

from . import __version__, appname


@clicommand(message="Calculating fitness with the following parameters")
def _fitness(args: Namespace) -> None:

    from .io import load_anndata
    from .models.anndata import AnnDataWLSModel
    from .transforms import compute_log_ratios

    print_err(
        f"Loading from {args.inputs}"
    )
    adata = load_anndata(
        counts=args.inputs,
        sample_meta=args.sample_sheet,
        strain_meta=args.barcode_sheet,
        reference=args.reference,
        count_column=args.count_column,
        timepoint_column=args.timepoint_column,
        t0=args.t0,
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
    model = AnnDataWLSModel()
    results = model.fit(adata=adata)
    print_err(results)
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
            pseudocount,
            plots,
            plot_format,
            formatting,
        ],
    )
    app = CLIApp(
        name=appname,
        version=__version__,
        description="Analysis of barcoded, pooled functional genomics assays.",
        commands=[
            fitness,
        ]
    )
    app.run()
    return None


if __name__ == '__main__':
    main()
