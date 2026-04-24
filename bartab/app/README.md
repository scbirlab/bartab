---
title: bartab
emoji: 🍹
colorFrom: blue
colorTo: green
sdk: gradio
sdk_version: "5.50.0"
app_file: app.py
pinned: true
short_description: Fitness from sequencing-based pooled competition experiments.
tags:
  - biology
  - sequencing
  - pooled-screen
  - fitness
  - gradio
---

# 🍹 bartab

**Estimate strain fitness from sequencing-based pooled competition experiments.**

`bartab` analyses pooled barcode sequencing experiments in which multiple strains,
mutants, guides, or constructs are grown together and quantified by NGS over time.
It estimates the relative fitness of each barcode compared with a reference strain,
typically WT.

[![Open in Spaces](https://huggingface.co/datasets/huggingface/badges/resolve/main/open-in-hf-spaces-md-dark.svg)](https://huggingface.co/spaces/scbirlab/bartab)

---

## What this app does

Upload three tables:

1. **Count table**  
   Read or UMI counts for each barcode in each sequencing sample.

2. **Sample sheet**  
   Metadata describing each sample: sample ID, culture/replicate ID, timepoint,
   optional concentration, optional growth measurement, and optional sampled volume.

3. **Barcode sheet**  
   Barcode or strain identifiers, plus optional strain-level metadata.

`bartab` then estimates how each barcode changes relative to a reference barcode as
the culture expands.

---

### Option 1: Spike-in normalisation

Use this when your experiment includes a non-growing spike-in barcode, such as a
heat-killed strain, plasmid-only control, or other fitness-zero control.

If the spike-in does not grow, then changes in the spike-in:WT count ratio report how
much WT has expanded.

For a non-growing spike-in:

$$
\log \left(
\frac{c_{spike}(t)}{c_{wt}(t)}
\frac{c_{wt}(0)}{c_{spike}(0)}
\right)
=
-
\log
\frac{n_{wt}(t)}{n_{wt}(0)}
$$

Substituting this into the barcode model gives:

$$
\log \left(
\frac{c_i(t)}{c_{wt}(t)}
\frac{c_{wt}(0)}{c_i(0)}
\right)
=
\left(
1 - \frac{w_i}{w_{wt}}
\right)
\log \left(
\frac{c_{spike}(t)}{c_{wt}(t)}
\frac{c_{wt}(0)}{c_{spike}(0)}
\right)
$$

So, when using spike-in normalisation, `bartab` estimates relative fitness from the
relationship between:

- the barcode:WT log-ratio change
- the spike-in:WT log-ratio change

---

### Option 2: Growth-based normalisation

Use this when you have measured culture growth directly, for example:

- OD600
- CFU/mL
- estimated generations
- another density-like growth measurement

In this mode, `bartab` uses the supplied growth column as the culture expansion axis
instead of estimating expansion from a spike-in barcode.

---

## Analysis modes

### 👟 Fitness in a single condition

This is the standard mode for pooled competition assays.

`bartab` fits a weighted least-squares model to estimate relative fitness for each
barcode.

Main outputs:

| Column | Meaning |
|---|---|
| `fitness` | relative fitness compared with the reference barcode |
| `fitness_low` | lower confidence interval bound |
| `fitness_high` | upper confidence interval bound |
| `slope_p` | p-value for deviation from neutral fitness |

Interpretation:

| Fitness | Interpretation |
|---|---|
| `1` | grows like the reference |
| `< 1` | growth disadvantage |
| `> 1` | growth advantage |

---

### 📉 Dose response

Use this when the sample sheet contains a concentration column.

`bartab` first estimates barcode fitness across concentrations, then fits a
two-parameter Hill model to estimate dose-response parameters.

Main outputs:

| Column | Meaning |
|---|---|
| `log_ic50` | log10 concentration giving 50% inhibition |
| `log_ic50_p` | p-value associated with the IC50 estimate |

Lower IC50 values indicate greater sensitivity to the inducer or drug.

---

## Required input tables

The app lets you select the relevant columns after upload.

Supported file formats:

- `.csv`
- `.tsv`
- `.txt`
- `.xlsx`

---

### 1. Count table

One row per barcode per sample.

Required columns:

| Meaning | Example |
|---|---|
| barcode / strain identifier | `strain_id` |
| sample identifier | `sample_id` |
| read or UMI count | `count` |

Example:

| strain_id | sample_id | count |
|---|---:|---:|
| wt | sample_0 | 12034 |
| mutant_A | sample_0 | 8312 |
| spike | sample_0 | 5021 |
| wt | sample_1 | 18420 |
| mutant_A | sample_1 | 6420 |
| spike | sample_1 | 2100 |

---

### 2. Barcode sheet

One row per barcode.

Required columns:

| Meaning | Example |
|---|---|
| barcode / strain identifier | `strain_id` |

Optional metadata columns are carried through to the output.

Example:

| strain_id | gene | annotation |
|---|---|---|
| wt | WT | reference |
| mutant_A | geneA | deletion mutant |
| spike | spike | non-growing spike-in |

---

### 3. Sample sheet

One row per sequencing sample.

Required columns:

| Meaning | Example |
|---|---|
| sample identifier | `sample_id` |
| culture / biological replicate | `replicate` |
| timepoint | `timepoint` |

Optional columns:

| Meaning | Example | Used for |
|---|---|---|
| concentration | `dose` | dose-response analysis |
| growth measurement | `growth` | growth-based normalisation |
| sampled volume | `volume` | adaptive-volume sampling |

Example:

| sample_id | replicate | timepoint | dose | growth | volume |
|---|---|---:|---:|---:|---:|
| sample_0 | rep1 | 0 | 0 | 0.05 | 1.0 |
| sample_1 | rep1 | 1 | 0 | 0.20 | 1.0 |
| sample_2 | rep1 | 2 | 0 | 0.80 | 1.0 |

---

## Example datasets

The Space includes synthetic example datasets for:

- single-concentration analysis using spike-in normalisation
- single-concentration analysis using growth-based normalisation
- dose-response analysis using spike-in normalisation
- dose-response analysis using growth-based normalisation

Start with these examples if you want to see the expected table structure.

---

## Outputs

After analysis, the app returns:

1. **Fitted parameter table**  
   A CSV file containing estimated fitness or dose-response parameters.

2. **Annotated `.h5ad` object**  
   A full `AnnData` object containing input data, metadata, transformations,
   fitted values, and model outputs.

3. **Diagnostic plots**

Depending on the analysis mode, plots include:

- time vs count
- expansion vs count
- time vs ratio
- expansion vs ratio
- predicted vs observed
- volcano plot
- dose-response curves

---

## Practical guidance

- Use consistent barcode and sample identifiers across all three tables.
- If using spike-in normalisation, include the spike-in barcode in both the count table
and barcode sheet.
- If using growth-based normalisation, include a numeric growth column in the sample sheet.
- For dose-response analysis, concentration values must be numeric.
- Low-count barcodes may give unstable estimates.
- Always inspect diagnostic plots before interpreting individual hits.

---

## The model

[Interactive tutorial on analysis principles](https://huggingface.co/spaces/scbirlab/tutorial-seq-fitness)

In a pooled growth experiment, strains compete in the same culture. Their absolute
growth curves may be complex because the pool eventually approaches carrying capacity.
However, if we compare every strain to a reference strain, the shared density-dependent
term cancels. The reference strain’s expansion can then be used as the effective
growth clock.

For strain $i$ relative to WT:

$$
\log n_i(t)
=
\frac{w_i}{w_{wt}}
\log \frac{n_{wt}(t)}{n_{wt}(0)}
+
\log n_i(0)
$$

where:

- $n_i(t)$ is the abundance of strain $i$ at time $t$
- $w_i$ is the intrinsic growth rate of strain $i$
- $w_i / w_{wt}$ is the relative fitness

So the problem becomes: estimate how barcode abundance changes relative to WT as WT
expands.

---

## From cells to sequencing counts

In practice, we do not observe cell numbers directly. We observe read counts:

$$
c_i(t)
$$

These are affected by sampling, library preparation, and sequencing depth. To remove
sample-specific sequencing depth effects, `bartab` works with ratios of barcode counts
to the reference barcode.

The key quantity is:

$$
\log \left(
\frac{c_i(t)}{c_{wt}(t)}
\frac{c_{wt}(0)}{c_i(0)}
\right)
$$

This is the log-change in the abundance of barcode $i$ relative to WT, normalised
to the starting timepoint.

Under the model:

$$
\log \left(
\frac{c_i(t)}{c_{wt}(t)}
\frac{c_{wt}(0)}{c_i(0)}
\right)
=
\left(
\frac{w_i}{w_{wt}} - 1
\right)
\log
\frac{n_{wt}(t)}{n_{wt}(0)}
$$

Thus, each barcode should follow an approximately straight line. The slope gives the
barcode’s fitness relative to WT.

---

## Estimating culture expansion

The remaining problem is that the true WT expansion,

$$
\frac{n_{wt}(t)}{n_{wt}(0)}
$$

is usually not directly observed. `bartab` supports two ways to estimate it.

---

## When this model is appropriate

`bartab` is designed for pooled competition experiments where:

- barcodes identify strains, mutants, guides, or constructs
- all barcodes are grown together in shared cultures
- barcode abundance is quantified by sequencing
- one barcode can be treated as a reference
- expansion can be estimated from a spike-in or measured growth

It is especially useful for:

- bacterial pooled competition assays
- barcoded mutant libraries
- CRISPRi or guide-based growth assays
- chemical-genetic pooled fitness experiments
- inducer or drug dose-response screens

---

## Limitations

`bartab` estimates relative fitness, not absolute growth rate.

Results may be unreliable when:

- the reference barcode is depleted or poorly counted
- the spike-in is not truly non-growing
- barcode counts are extremely low
- barcode identities are mismatched between tables
- bottlenecks dominate the experiment
- strong barcode-specific sequencing biases are present
- timepoints or concentrations are too sparse for model fitting

For dose-response fitting, IC50 estimates are most meaningful when the tested
concentration range brackets the transition from weak to strong inhibition.

---

## Local use

To run the app locally:

```bash
pip install -r requirements.txt
gradio app.py
```

For package and source code, see https://github.com/scbirlab/bartab
