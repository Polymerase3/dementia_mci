# Dementia PhIP-Seq analysis pipeline

[![phiper](https://img.shields.io/badge/phiper-Polymerase3%2Fphiper-blue)](https://github.com/Polymerase3/phiper)
[![Lifecycle: stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: inactive](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

This repository contains a small, frozen analysis pipeline for PhIP-Seq dementia data.
The goal is reproducibility, not ongoing development.

---

## 1. What you need

- A working installation of **R** (≥ 4.x).
- Command-line access (Linux/macOS terminal or Windows PowerShell).
- This repository cloned locally.
- The `phiper` source tarball is already included in the repository root (for example: `phiper_x.y.z.tar.gz`).

You do **not** need to download any external data: all required input files are part of the repo or created by the scripts.

---

## 2. Clone the repository

Open a terminal and run:

```bash
git clone https://github.com/Polymerase3/dementia_mci.git
cd dementia_mci
```

All remaining commands below assume you are inside this project directory.

---

## 3. Install `phiper` from the bundled tarball

From the project root (the folder where you see the `phiper_*.tar.gz` file), install `phiper` directly from the tarball.

### Option A: Command-line install (recommended)

```bash
R CMD INSTALL phiper_*.tar.gz
```

The `*` will expand to the actual versioned tarball name (for example `phiper_0.3.0.tar.gz`).

### Option B: From within an R session

If you prefer to work inside R:

```r
install.packages("phiper_0.3.0.tar.gz", repos = NULL, type = "source")
```

Make sure you run this **from the project root**, or provide the full path to the tarball.

Once this completes without errors, the `phiper` package is available for the scripts.

---

## 4. Run the pipeline

There are two scripts, and you should run them **in order**.

### 4.1. Step 1 – Data preparation

This script prepares the input data (for example, the Parquet file used by the analysis):

```bash
Rscript 01-data_prepare.R
```

- Run this once.
- It will create the intermediate data files expected by the main analysis script
  (e.g. `data_original/other_data/dementia_full.parquet`).

### 4.2. Step 2 – Main analysis

The main script performs all analyses (alpha/beta diversity, POP framework, DELTA framework, and all plots).
It writes everything into a newly created `results/` directory.

**Important:** The script expects arguments in the form `KEY=VALUE`. You should call it like this:

```bash
Rscript R/02-analysis_script.R N_CORES=30 MAX_GB=50 LOG=TRUE
```

#### Arguments

- `N_CORES`  
  Number of parallel workers used by the DELTA framework and some other steps.  
  Example: `N_CORES=30`  
  Choose this according to how many CPU cores you actually have or are allowed to use.

- `MAX_GB`  
  Rough memory limit (in gigabytes) for `future.globals.maxSize`.  
  Default in the script is `40`.  
  Example: `MAX_GB=50`.

- `LOG`  
  Controls logging for the DELTA framework step (`ph_prevalence_shift2`).  
  - `LOG=TRUE`  – write a log file per comparison.
  - `LOG=FALSE` – no log files.

- `LOG_FILE` (optional)  
  Path to a single log file that should be used instead of the per-comparison default.  
  If you do not provide this, each comparison gets its own  
  `results/<group1>_vs_<group2>/DELTA_framework/log.txt`.

Example with all parameters:

```bash
Rscript R/02-analysis_script.R N_CORES=24 MAX_GB=60 LOG=TRUE LOG_FILE="my_delta_log.txt"
```

If everything runs correctly, you will end up with a populated `results/` directory as described below.

---

## 5. Output: structure of the `results/` directory

The main script creates a top-level `results/` folder and then one subfolder per comparison.

### 5.1. Top level

```text
results/
  peptide_library.rds
  <group1>_vs_<group2>/
  <group1>_vs_<group2>/
  ...
```

- `peptide_library.rds`  
  RDS with the full peptide library as a plain `data.frame`.  
  This is used by several downstream steps (POP and DELTA frameworks).

Each comparison (e.g. `control_vs_dementia`, `control_vs_MCI`) has its own directory:

```text
results/
  control_vs_dementia/
  control_vs_MCI/
  ...
```

### 5.2. Per-comparison directory

Inside each comparison directory (here shown for `control_vs_dementia/`):

```text
results/
  control_vs_dementia/
    control_vs_dementia_data.rds
    enrichment_counts.svg
    alpha_diversity/
    beta_diversity/
    POP_framework/
    DELTA_framework/
```

#### Root files

- `control_vs_dementia_data.rds`  
  Filtered PHIP-Seq data for this comparison (only the two groups involved), saved as a `data.frame`.

- `enrichment_counts.svg`  
  Static SVG plot from `plot_enrichment_counts()`: enrichment counts by group.

---

### 5.3. `alpha_diversity/`

```text
results/<cmp>/alpha_diversity/
  table.xlsx
  plot.svg
```

- `table.xlsx`  
  Excel table with alpha diversity metrics (e.g. richness) per sample and group.

- `plot.svg`  
  Static alpha diversity plot from `plot_alpha_diversity()`, showing differences between groups.

---

### 5.4. `beta_diversity/`

```text
results/<cmp>/beta_diversity/
  distance_matrix.xlsx
  pcoa_results.rds
  capscale_results.rds
  permanova_results.rds
  dispersion_results.rds
  tsne2d_results.xlsx
  tsne2d_plot.svg
  tsne3d_results.xlsx
  tsne3d_plot.html
  pcoa_plot.svg
  scree_plot.svg
  dispersion_plot.svg
```

- `distance_matrix.xlsx`  
  Bray–Curtis distance matrix (after Hellinger normalization) as an Excel file.

- `pcoa_results.rds`  
  R object from `compute_pcoa()` with eigenvalues, sample coordinates, and feature loadings.

- `capscale_results.rds`  
  R object from `compute_capscale()` (constrained ordination / dbRDA for `group_char`).

- `permanova_results.rds`  
  R object from `compute_permanova()` with PERMANOVA results and pairwise contrasts.

- `dispersion_results.rds`  
  R object from `compute_dispersion()` with group dispersions and distances.

- `tsne2d_results.xlsx`  
  t-SNE coordinates in 2D with metadata.

- `tsne2d_plot.svg`  
  Static 2D t-SNE plot (`plot_tsne(view = "2d")`).

- `tsne3d_results.xlsx`  
  t-SNE coordinates in 3D with metadata.

- `tsne3d_plot.html`  
  Interactive 3D t-SNE widget (open in a web browser).

- `pcoa_plot.svg`  
  PCoA plot with group centroids and ellipses.

- `scree_plot.svg`  
  Scree plot for the first 15 PCoA axes.

- `dispersion_plot.svg`  
  Plot of group dispersions (violin + box + points) for the selected contrast.

---

### 5.5. `POP_framework/`

```text
results/<cmp>/POP_framework/
  single_peptide.csv
  taxa_ranks.csv
  plots/
    peptide_id_static.svg
    peptide_id_interactive.html
    phylum_static.svg
    phylum_interactive.html
    class_static.svg
    class_interactive.html
    ...
```

- `single_peptide.csv`  
  Output of `ph_prevalence_compare()` at the peptide level (`rank_cols = "peptide_id"`).

- `taxa_ranks.csv`  
  Output of `ph_prevalence_compare()` aggregated at different taxonomic ranks  
  (`phylum`, `class`, `order`, `family`, `genus`, `species`).

- `plots/`  
  For each rank (including `peptide_id`) there are:

  - `<rank>_static.svg` – static scatter plot (`scatter_static()`),  
    showing prevalence changes between the two groups.
  - `<rank>_interactive.html` – interactive version (`scatter_interactive()`),  
    openable in any web browser.

---

### 5.6. `DELTA_framework/`

```text
results/<cmp>/DELTA_framework/
  log.txt
  delta_table.csv
  uncorrected_significant_static_forestplot.svg
  BHcorrected_significant_static_forestplot.svg
  uncorrected_significant_interactive_forestplot.html
  BHcorrected_significant_interactive_forestplot.html
  interesting_features/
    <feature>_scatter_static.svg
    <feature>_scatter_interactive.html
    <feature>_deltaplot_static.svg
    <feature>_deltaplot_interactive.html
    <feature>_ecdfplot_static.svg
    <feature>_ecdfplot_interactive.html
```

- `log.txt`  
  Log file produced by `ph_prevalence_shift2()` if `LOG=TRUE` and `LOG_FILE` is not set.

- `delta_table.csv`  
  Main DELTA framework results table from `ph_prevalence_shift2()`, with:
  - test statistics (`T_obs`),
  - permutation p-values (`p_perm`),
  - adjusted p-values (`p_adj_rank`),
  - and an additional column `Z` (p-to-Z transform).

- `uncorrected_significant_static_forestplot.svg`  
  Static forest plot of features significant by raw `p_perm`, using `stat = "Z"`.

- `BHcorrected_significant_static_forestplot.svg`  
  Static forest plot of features significant after BH correction (`p_adj_rank`,  
  `sig_level = 0.10`), again using Z-scores.

- `uncorrected_significant_interactive_forestplot.html`  
  Interactive forest plot (raw `p_perm`), open in a browser.

- `BHcorrected_significant_interactive_forestplot.html`  
  Interactive forest plot (BH-adjusted, `p_adj_rank`), open in a browser.

#### `interesting_features/`

For all features with `p_perm < 0.05` that can be mapped to the peptide library,
the script creates a set of diagnostic plots:

- `<feature>_scatter_static.svg`  
  Static prevalence scatter for the feature (per-peptide prevalence in both groups).

- `<feature>_scatter_interactive.html`  
  Interactive version of the prevalence scatter.

- `<feature>_deltaplot_static.svg`  
  Static DELTA plot (`deltaplot_prevalence()`), with optional smoothing depending on
  the number of points.

- `<feature>_deltaplot_interactive.html`  
  Interactive DELTA plot (`deltaplot_prevalence_interactive()`).

- `<feature>_ecdfplot_static.svg`  
  Static ECDF comparison (`ecdf_prevalence()`), including medians and KS statistics.

- `<feature>_ecdfplot_interactive.html`  
  Interactive ECDF plot (`ecdf_prevalence_interactive()`).

---

## 6. Re-running and cleaning up

- Re-running `R/02-analysis_script.R` will **overwrite** the contents of the
  corresponding comparison directories inside `results/` if they already exist.
- If you want to keep old outputs, copy the whole `results/` folder elsewhere
  before running the pipeline again.
