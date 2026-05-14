# 🍹 bartab

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/bartab/python-test.yml)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/bartab)
![PyPI](https://img.shields.io/pypi/v/bartab)
[![Open in Spaces](https://huggingface.co/datasets/huggingface/badges/resolve/main/open-in-hf-spaces-md-dark.svg)](https://huggingface.co/spaces/scbirlab/bartab)

<img src="docs/source/_static/logo.svg" width="200">

**bartab** estimates relative fitness from sequencing-based pooled competition experiments.

It is designed for experiments where barcoded strains, guides, mutants, or constructs are grown together, sampled over time, and quantified by NGS read counts or UMI counts. `bartab` estimates the fitness of each barcode relative to a reference barcode, usually WT, using either spike-in normalisation or measured culture growth.

`bartab` provides:

- a command-line interface for fitting, plotting, and simulating pooled fitness experiments
- a Python API built around `AnnData`
- weighted least-squares fitness estimation
- dose-response fitting with a Hill/log-logistic model
- diagnostic plots
- synthetic data simulation

---

## Contents

- [Installation](#installation)
- [Conceptual model](#conceptual-model)
- [Input tables](#input-tables)
- [Command-line interface](#command-line-interface)
  - [`bartab fit`](#bartab-fit)
  - [`bartab plot`](#bartab-plot)
  - [`bartab sim`](#bartab-sim)
- [Python API](#python-api)
  - [Load data](#load-data)
  - [Compute log-ratios](#compute-log-ratios)
  - [Fit fitness models](#fit-fitness-models)
  - [Fit dose-response models](#fit-dose-response-models)
  - [Plot results](#plot-results)
  - [Simulate data](#simulate-data)
- [Output structure](#output-structure)
- [Interpretation](#interpretation)
- [Limitations](#limitations)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Documentation](#documentation)

---

## Installation

### From PyPI

```bash
pip install bartab
```

### From source

```bash
git clone https://github.com/scbirlab/bartab.git
cd bartab
pip install -e .
```

For development, install test and documentation dependencies as appropriate for your local setup.

---

## Conceptual model

[Interactive tutorial on analysis principles](https://huggingface.co/spaces/scbirlab/tutorial-seq-fitness).

---

## Input tables

`bartab` expects three tables:

1. count table
2. sample sheet
3. barcode sheet

Files may be CSV, TSV, TXT, or other formats supported by `carabiner.pd.read_table`.

---

### Count table

One row per barcode per sample.

Required columns:

| Meaning | Default CLI option | Typical name |
|---|---:|---|
| barcode / strain identifier | `--barcode-column` | `strain_id` |
| sample identifier | `--sample-column` | `sample_id` |
| count value | `--count-column` | `count` |

Example:

| strain_id | sample_id | count |
|---|---:|---:|
| wt | sample_0 | 12034 |
| mutant_A | sample_0 | 8312 |
| spike | sample_0 | 5021 |
| wt | sample_1 | 18420 |
| mutant_A | sample_1 | 6420 |
| spike | sample_1 | 2100 |

Counts must be non-negative.

---

### Sample sheet

One row per sequencing sample.

Required columns:

| Meaning | Default CLI option | Typical name |
|---|---:|---|
| sample identifier | `--sample-column` | `sample_id` |
| culture / biological replicate | `--culture-column` | `culture_id` |
| timepoint | `--timepoint-column` | `timepoint` |

Optional columns:

| Meaning | CLI option | Used for |
|---|---:|---|
| concentration | `--concentration-column` | dose-response analysis |
| growth measurement | `--growth-column` | growth-based normalisation |
| sampled volume | `--volume-column` | adaptive-volume correction with spike-in |

Example:

| sample_id | culture_id | timepoint | dose | growth | volume |
|---|---|---:|---:|---:|---:|
| sample_0 | rep1 | 0 | 0 | 0.05 | 1.0 |
| sample_1 | rep1 | 1 | 0 | 0.20 | 1.0 |
| sample_2 | rep1 | 2 | 0 | 0.80 | 1.0 |

If `--t0` is not supplied, the minimum value in the timepoint column is used as the initial timepoint.

---

### Barcode sheet

One row per barcode.

Required columns:

| Meaning | Default CLI option | Typical name |
|---|---:|---|
| barcode / strain identifier | `--barcode-column` | `strain_id` |

Optional metadata columns are preserved in the output `AnnData.obs`.

Example:

| strain_id | gene | annotation |
|---|---|---|
| wt | WT | reference |
| mutant_A | geneA | deletion mutant |
| spike | spike | non-growing spike-in |

---

## Command-line interface

Run:

```bash
bartab --help
```

`bartab` has three main commands:

```bash
bartab fit   # estimate fitness
bartab plot  # plot fitted results
bartab sim   # simulate pooled competition data
```

---

## `bartab fit`

Estimate relative fitness per strain or barcode.

### Minimal spike-in example

```bash
bartab fit counts.csv \
  --sample-sheet sample_meta.csv \
  --barcode-sheet strain_meta.csv \
  --output results.h5ad \
  --reference wt \
  --spike-name spike \
  --use-spike \
  --barcode-column strain_id \
  --sample-column sample_id \
  --culture-column culture_id \
  --count-column count \
  --timepoint-column timepoint
```

This writes:

```text
results.h5ad
results.csv
```

when the output path ends in `.h5ad`.

---

### Minimal growth-based example

```bash
bartab fit counts.csv \
  --sample-sheet sample_meta.csv \
  --barcode-sheet strain_meta.csv \
  --output results.h5ad \
  --reference wt \
  --barcode-column strain_id \
  --sample-column sample_id \
  --culture-column culture_id \
  --count-column count \
  --timepoint-column timepoint \
  --growth-column growth \
  --growth-type density
```

---

### Dose-response example

```bash
bartab fit counts.csv \
  --sample-sheet sample_meta.csv \
  --barcode-sheet strain_meta.csv \
  --output dose_response.h5ad \
  --reference wt \
  --spike-name spike \
  --use-spike \
  --barcode-column strain_id \
  --sample-column sample_id \
  --culture-column culture_id \
  --count-column count \
  --timepoint-column timepoint \
  --concentration-column dose
```

If `--concentration-column` is supplied, `bartab` first fits per-concentration WLS fitness estimates and then fits a Hill/log-logistic dose-response model.

---

### Output format

If `--output` ends in:

| Output suffix | Behaviour |
|---|---|
| `.h5ad` | writes full annotated `AnnData` object and a companion `.csv` table |
| `.h5ad.gz` | writes compressed `AnnData` and companion `.csv` table |
| `.csv` | writes fitted parameter table only |
| `.tsv` / `.txt` | writes fitted parameter table only, tab-separated |

---

### Main `fit` options

| Option | Default | Description |
|---|---:|---|
| `input` | required | count table |
| `--output`, `-o` | required | output file |
| `--sample-sheet`, `-c` | required | sample metadata table |
| `--barcode-sheet`, `-b` | required | barcode metadata table |
| `--reference`, `-r` | required | reference barcode/strain |
| `--barcode-column`, `-u` | `strain_name` | barcode column in count and barcode tables |
| `--sample-column`, `-s` | `sample_id` | sample ID column |
| `--culture-column` | `culture_id` | biological replicate/culture column |
| `--count-column`, `-a` | `count` | count column |
| `--timepoint-column`, `-t` | `timepoint` | timepoint column |
| `--t0`, `-0` | minimum timepoint | initial timepoint value |
| `--spike-name`, `-x` | `None` | spike-in barcode name |
| `--use-spike` | `False` | use spike-in to estimate culture expansion |
| `--growth-column`, `-d` | `OD600` | growth measurement column |
| `--growth-type` | `density` | either `density` or `generations` |
| `--volume-column` | `None` | sampled volume column |
| `--concentration-column`, `-k` | `None` | concentration column for dose-response |
| `--pseudocount` | `1.0` | value added to counts before log transform |
| `--model-type` | `WLS` | currently accepted: `WLS`, `OLS`, `HillFitnessModel` |

---

## `bartab plot`

Generate diagnostic plots from a fitted `.h5ad` file.

```bash
bartab plot results.h5ad \
  --output results_plots \
  --model-type WLS \
  --plot-format png
```

This writes plots using the supplied prefix:

```text
results_plots_time-count.png
results_plots_time-ratio.png
results_plots_expansion-ratio.png
results_plots_pred-obs.png
results_plots_volcano.png
```

For dose-response/Hill results:

```bash
bartab plot dose_response.h5ad \
  --output dose_response_plots \
  --model-type HillFitnessModel \
  --plot-format png
```

This additionally writes:

```text
dose_response_plots_dr.png
```

---

### Plot options

| Option | Default | Description |
|---|---:|---|
| `input` | required | fitted `.h5ad` file |
| `--output`, `-o` | required | filename prefix |
| `--highlight`, `-X` | `None` | barcode names to highlight |
| `--model-type` | `WLS` | model to plot; `WLS`, `OLS`, or `HillFitnessModel` |
| `--plot-format` | `png` | `png`, `PNG`, `pdf`, or `PDF` |

---

## `bartab sim`

Simulate count data from a pooled competition experiment.

The simulator takes a JSON file describing true strain fitness values and writes three CSV files:

```text
<output>_count.csv
<output>_sample_meta.csv
<output>_strain_meta.csv
```

---

### Single-concentration simulation

Input JSON:

```json
{
  "wt": 1.0,
  "spike": 0.0,
  "mutant_fast": 1.25,
  "mutant_slow": 0.5
}
```

Run:

```bash
bartab sim fitness.json \
  --output synthetic/single \
  --timepoints 8 \
  --n-cultures 3 \
  --reads-per-barcode 1000 \
  --seed 42
```

Then fit:

```bash
bartab fit synthetic/single_count.csv \
  --sample-sheet synthetic/single_sample_meta.csv \
  --barcode-sheet synthetic/single_strain_meta.csv \
  --output synthetic/single_results.h5ad \
  --reference wt \
  --spike-name spike \
  --use-spike \
  --barcode-column strain_id \
  --sample-column sample_id \
  --culture-column replicate \
  --count-column count \
  --timepoint-column timepoint
```

---

### Dose-response simulation

For dose-response simulation, the JSON values are interpreted as `[IC50, bottom]`.

Input JSON:

```json
{
  "wt": [10000.0, 1.0],
  "spike": [10000.0, 0.0],
  "sensitive": [10.0, 0.0],
  "resistant": [1000.0, 0.0]
}
```

Run:

```bash
bartab sim dose_fitness.json \
  --output synthetic/dose \
  --n-dose 6 \
  --dose-max 1000 \
  --dose-fold 2 \
  --timepoints 8 \
  --n-cultures 3 \
  --reads-per-barcode 1000 \
  --seed 42
```

Then fit:

```bash
bartab fit synthetic/dose_count.csv \
  --sample-sheet synthetic/dose_sample_meta.csv \
  --barcode-sheet synthetic/dose_strain_meta.csv \
  --output synthetic/dose_results.h5ad \
  --reference wt \
  --spike-name spike \
  --use-spike \
  --barcode-column strain_id \
  --sample-column sample_id \
  --culture-column replicate \
  --count-column count \
  --timepoint-column timepoint \
  --concentration-column dose
```

---

### Simulation options

| Option | Default | Description |
|---|---:|---|
| `input` | required | JSON file of strain fitness values |
| `--output`, `-o` | required | output prefix |
| `--generate-controls`, `-z` | `0` | number of neutral controls to generate |
| `--generate-more`, `-m` | `0` | number of additional random-fitness strains |
| `--inoculum`, `-n` | `1000` | cells per strain in inoculum |
| `--carrying-capacity`, `-K` | `10.0` | maximum fold-change supported |
| `--timepoints`, `-s` | `10` | number of sampled timepoints |
| `--max-time`, `-t` | `10` | final simulated timepoint |
| `--n-cultures`, `-r` | `3` | biological replicate cultures |
| `--reads-per-barcode`, `-b` | `1000` | mean sequencing depth per barcode per sample |
| `--seed`, `-e` | `42` | random seed |
| `--n-dose`, `-d` | `None` | number of dose-response concentrations |
| `--dose-max` | `1000` | maximum simulated dose |
| `--dose-fold` | `2` | dilution spacing parameter |

---

## Python API

The Python API is built around `AnnData`.

A typical workflow is:

1. load count, sample, and barcode tables into `AnnData`
2. compute barcode/reference log-ratios and expansion axis
3. fit a model
4. inspect `adata.obs`
5. optionally plot or write output

---

## Load data

```python
from bartab.io import load_anndata

adata = load_anndata(
    counts="counts.csv",
    sample_meta="sample_meta.csv",
    strain_meta="strain_meta.csv",
    reference="wt",
    spike="spike",
    count_column="count",
    strain_id="strain_id",
    sample_id="sample_id",
    culture_id="culture_id",
    timepoint_column="timepoint",
)
```

The same function also accepts pandas `DataFrame` objects:

```python
adata = load_anndata(
    counts=counts_df,
    sample_meta=sample_df,
    strain_meta=strain_df,
    reference="wt",
    spike="spike",
    count_column="count",
    strain_id="strain_id",
    sample_id="sample_id",
    culture_id="culture_id",
    timepoint_column="timepoint",
)
```

For dose-response data:

```python
adata = load_anndata(
    counts=counts_df,
    sample_meta=sample_df,
    strain_meta=strain_df,
    reference="wt",
    spike="spike",
    count_column="count",
    strain_id="strain_id",
    sample_id="sample_id",
    culture_id="culture_id",
    timepoint_column="timepoint",
    concentration_column="dose",
)
```

---

## Compute log-ratios

### Spike-in normalisation

```python
from bartab.transforms import compute_log_ratios

adata = compute_log_ratios(
    adata,
    pseudocount=1.0,
    use_spike=True,
)
```

With adaptive-volume correction:

```python
adata = compute_log_ratios(
    adata,
    pseudocount=1.0,
    use_spike=True,
    volume_column="volume",
)
```

### Growth-based normalisation

```python
adata = compute_log_ratios(
    adata,
    pseudocount=1.0,
    growth_column="growth",
    growth_type="density",
    use_spike=False,
)
```

If the growth column already contains generations:

```python
adata = compute_log_ratios(
    adata,
    pseudocount=1.0,
    growth_column="generations",
    growth_type="generations",
    use_spike=False,
)
```

This adds:

```python
adata.layers["__log_ratio__"]
adata.var["__log_expansion__"]
```

---

## Fit fitness models

### Weighted least squares

```python
from bartab.models.anndata import AnnDataWLSModel

model = AnnDataWLSModel()
adata = model.fit(adata)
```

Results are written into `adata.obs` with names prefixed by the model name.

For WLS, common output columns include:

```text
WLS:slope
WLS:slope_p
WLS:slope_se
WLS:slope_ci_low
WLS:slope_ci_high
WLS:fitness
WLS:fitness_low
WLS:fitness_high
WLS:nobs
WLS:rsq
WLS:fit_status
```

Predictions are written to:

```python
adata.layers["WLS:predicted"]
```

The model name is also appended to:

```python
adata.uns["models_fitted"]
```

---

### Ordinary least squares

```python
from bartab.models.anndata import AnnDataOLSModel

model = AnnDataOLSModel()
adata = model.fit(adata)
```

OLS is mainly useful as an unweighted baseline or diagnostic comparison.

---

## Fit dose-response models

For dose-response analysis, first load data with `concentration_column` set, then compute log-ratios, fit WLS, and fit the Hill model.

```python
from bartab.io import load_anndata
from bartab.transforms import compute_log_ratios
from bartab.models.anndata import AnnDataWLSModel, AnnDataHillModel

adata = load_anndata(
    counts="dose_count.csv",
    sample_meta="dose_sample_meta.csv",
    strain_meta="dose_strain_meta.csv",
    reference="wt",
    spike="spike",
    count_column="count",
    strain_id="strain_id",
    sample_id="sample_id",
    culture_id="replicate",
    timepoint_column="timepoint",
    concentration_column="dose",
)

adata = compute_log_ratios(
    adata,
    pseudocount=1.0,
    use_spike=True,
)

adata = AnnDataWLSModel().fit(adata)

adata = AnnDataHillModel().fit(
    adata,
    concentration="dose",
)
```

Hill model outputs include:

```text
HillFitnessModel:log_ic50
HillFitnessModel:log_ic50_p
HillFitnessModel:log_ic50_se
HillFitnessModel:log_ic50_ci_low
HillFitnessModel:log_ic50_ci_high
HillFitnessModel:ic50
HillFitnessModel:ic50_low
HillFitnessModel:ic50_high
HillFitnessModel:h
HillFitnessModel:h_p
HillFitnessModel:nobs
HillFitnessModel:dof
HillFitnessModel:fit_status
```

The Hill model uses a log-logistic form fit on log concentration. Concentration-zero samples are omitted from the non-linear fit.

---

## Plot results

Plotting functions live in `bartab.plotting`.

```python
from bartab.plotting import (
    time_vs_count,
    time_vs_ratio,
    expansion_vs_count,
    expansion_vs_ratio,
    pred_vs_true,
    volcano,
    dose_response,
)
```

Examples:

```python
fig, ax = time_vs_count(
    adata,
    highlight_barcodes=["mutant_A"],
    filename="plots/time_count.png",
)
```

```python
fig, ax = expansion_vs_ratio(
    adata,
    highlight_barcodes=["mutant_A"],
    filename="plots/expansion_ratio.png",
)
```

```python
fig, ax = pred_vs_true(
    adata,
    model_name="WLS",
    highlight_barcodes=["mutant_A"],
    filename="plots/pred_obs.png",
)
```

```python
fig, ax = volcano(
    adata,
    model_name="WLS",
    param="fitness",
    p="slope_p",
    vline=1.0,
    highlight_barcodes=["mutant_A"],
    filename="plots/volcano.png",
)
```

For dose-response:

```python
fig, ax = dose_response(
    adata,
    model_name="WLS",
    highlight_barcodes=["sensitive", "resistant"],
    filename="plots/dose_response.png",
)
```

Plot files are saved at high resolution. When a `filename` is provided, `bartab` also writes the plotted data alongside the figure through the internal `figsaver` utility.

---

## Simulate data

The simulator can be used from Python.

```python
from bartab.simulation import calculate_growth_curves, reads_sampler

fitness = {
    "wt": 1.0,
    "mutant_slow": 0.5,
    "mutant_fast": 1.25,
}

t, ref_expansion, growths = calculate_growth_curves(
    inoculum=1_000,
    fitness=fitness,
    inoculum_var=0.1,
    carrying_capacity=10,
    n_timepoints=8,
    max_time=10.0,
    ref_key="wt",
    seed=42,
)

counts = reads_sampler(
    growths,
    seq_depth=1_000,
    sample_frac=0.1,
    reps=3,
    variance=0.001,
    seed=42,
)
```

`counts` has shape:

```text
n_strains × n_replicates × n_timepoints
```

---

## Low-level model API

The lower-level model classes operate directly on arrays.

```python
import numpy as np
from bartab.models.linear import OLSModel, WLSModel

Y = np.array([
    [0.0, -0.5, -1.0],
    [0.0,  0.0,  0.0],
])

x = np.array([0.0, 1.0, 2.0])

model = OLSModel()
results, preds = model.fit(Y, x)
```

`results` is a list of dictionaries, one per row of `Y`.

For most users, the `AnnData` model API is preferable.

---

## Output structure

`bartab` stores results in an `AnnData` object.

### Main matrix

```python
adata.X
```

Raw count matrix, with:

```text
rows    = barcodes / strains
columns = samples
```

---

### Barcode metadata

```python
adata.obs
```

Contains barcode metadata and model outputs.

Important internal columns include:

| Column | Meaning |
|---|---|
| `__is_reference__` | whether the barcode is the reference |
| `__is_spike__` | whether the barcode is the spike-in |
| model-prefixed columns | fitted parameters and statistics |

Example model columns:

```text
WLS:fitness
WLS:fitness_low
WLS:fitness_high
WLS:slope_p
HillFitnessModel:ic50
HillFitnessModel:log_ic50_p
```

---

### Sample metadata

```python
adata.var
```

Contains sample metadata.

Important internal columns include:

| Column | Meaning |
|---|---|
| `__is_t0__` | whether the sample is an initial timepoint |
| `__culture_index__` | culture/replicate identifier |
| `__inducer__` | concentration or `"single dose"` |
| `__log_expansion__` | fitted x-axis: culture expansion |

---

### Layers

```python
adata.layers
```

Important layers include:

| Layer | Meaning |
|---|---|
| `__log_ratio__` | observed barcode/reference log-ratio change |
| `WLS:predicted` | WLS fitted values |
| `OLS:predicted` | OLS fitted values |
| `HillFitnessModel:predicted` | Hill model fitted values |

---

### Unstructured metadata

```python
adata.uns
```

Includes:

| Key | Meaning |
|---|---|
| `reference` | reference barcode name |
| `spike` | spike-in barcode name |
| `timepoint_column` | timepoint column used |
| `count_column` | count column used |
| `concentration_column` | concentration column, if any |
| `strain_id` | barcode ID column(s) |
| `sample_id` | sample ID column(s) |
| `culture_id` | culture/replicate column(s) |
| `models_fitted` | list of fitted models |

---

## Interpretation

For single-concentration WLS output:

| Output | Interpretation |
|---|---|
| `WLS:fitness = 1` | barcode grows like the reference |
| `WLS:fitness < 1` | growth disadvantage |
| `WLS:fitness > 1` | growth advantage |
| `WLS:slope_p` | evidence that the barcode deviates from neutral fitness |
| `WLS:fitness_low`, `WLS:fitness_high` | confidence interval for relative fitness |

For dose-response output:

| Output | Interpretation |
|---|---|
| `HillFitnessModel:ic50` | concentration giving 50% inhibition |
| lower IC50 | greater sensitivity |
| higher IC50 | greater resistance |
| `HillFitnessModel:h` | fitted Hill/logistic slope |
| `HillFitnessModel:log_ic50_p` | evidence for non-zero IC50 parameter estimate |

Dose-response estimates are most reliable when the tested concentrations bracket the fitness transition.

---

## Practical guidance

- Use the same barcode identifiers across the count table and barcode sheet.
- Use the same sample identifiers across the count table and sample sheet.
- If using spike-in normalisation, include exactly one spike-in barcode.
- If using growth-based normalisation, growth values must be positive.
- For dose-response analysis, concentrations must be numeric.
- Low-count barcodes can produce unstable estimates.
- Inspect diagnostic plots before interpreting individual hits.
- Check reference and spike-in behaviour before trusting global results.

---

## Limitations

Results may be unreliable when:

- the reference barcode is depleted or poorly counted
- the spike-in grows, dies, degrades, or is sampled inconsistently
- bottlenecks dominate the experiment
- counts are extremely low
- barcode identities are mismatched between tables
- dose-response concentrations do not span the active range

The WLS weights are estimated using an approximate delta-method variance on log 
count ratios with a method-of-moments dispersion estimate.

---

## Issues, problems, suggestions

Please add bug reports, feature requests, or suggestions to the issue tracker:

```text
https://www.github.com/scbirlab/bartab/issues
```

---

## Related resources

- Interactive Hugging Face Space: https://huggingface.co/spaces/scbirlab/bartab
- Tutorial Space: https://huggingface.co/spaces/scbirlab/tutorial-seq-fitness
- Source code: https://github.com/scbirlab/bartab
