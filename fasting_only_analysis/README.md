# Fasting-only bile acid analysis

This folder contains a separate longitudinal mixed-model analysis for fasting
conditions only:

- Baseline fasting: timepoint `2`
- NS4 fasting: timepoint `18`
- RTDS fasting: timepoint `26`

Timepoint `23` is excluded before model fitting. The analysis keeps the
corrected comparison workflow's core data handling:

- LTR median plate normalization for the 15 measured bile acids.
- Aggregate bile-acid rows remain missing when all component bile acids are
  missing.
- Positive adjusted concentrations are modeled as
  `log(concentration) ~ timepoint + (1 | subject)`.
- Satterthwaite degrees of freedom are used through `lmerTest` and `emmeans`.

## Why the FDR is recalculated

The fasting-only analysis answers a different hypothesis family from the
post-meal analysis. Excluding timepoint `23` refits the model using a different
set of observations and changes the multiple-testing set. The script therefore
recalculates raw p-values and Benjamini-Hochberg FDR q-values for fasting-only
tests rather than reusing q-values from the post-meal report.

The script reports two fasting-only result families:

1. Pairwise fasting contrasts:
   `18 vs 2`, `26 vs 2`, and `26 vs 18`.
2. A timepoint omnibus test asking whether any fasting timepoint differs.

For each result family, BH FDR is applied separately to:

- Individual bile acids.
- Aggregate classes.

## Run

From the repository root:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' fasting_only_analysis\fasting_only_lmm_analysis.R
```

## Outputs

The script writes generated outputs to `fasting_only_analysis/outputs`:

- `fasting_pairwise_results.csv`
- `fasting_omnibus_results.csv`
- `fasting_emm_results.csv`
- `verification_checks.csv`
- `session_info.txt`
- detailed fasting trajectory plots under `plots`
- grid plots for individual bile acids and aggregate classes under `plots`

## PCA

Run the fasting-only PCA from the repository root:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' fasting_only_analysis\fasting_only_pca.R
```

### Methodology

The PCA uses the same LTR median-adjusted fasting data, excludes timepoint
`23`, and uses `prcomp` on the natural-log transformed and standardized 15
measured bile acids. Aggregate classes are not used as PCA variables because
they are derived from the individual bile acids. The PCA is complete-case
across those 15 variables and does not impute missing values.

- **Variables**: 15 measured bile acids (aggregate classes excluded)
- **Preprocessing**: LTR median plate adjustment → natural-log → center + unit-variance scaling
- **Algorithm**: `prcomp(center = TRUE, scale. = TRUE)` via SVD
- **Missing values**: complete-case analysis (no imputation)

### PCA Outputs

The PCA script writes:

**PDF report** (`Fasting_Only_PCA_Report.pdf` in `outputs/`)

- Methods summary page
- Score plot with 68%/95% confidence ellipses and group centroids
- Subject trajectory plot (individual movement across timepoints)
- Scree plot + loading directions (side by side)
- Score distributions by condition (box + strip plots)
- Top signed loadings for PC1–PC3
- Variable contribution percentages
- Output file listing

**CSV files** (in `outputs/pca/`):

| File | Description |
|------|-------------|
| `fasting_pca_scores.csv` | PC scores for each retained observation |
| `fasting_pca_loadings.csv` | Variable loadings for all PCs |
| `fasting_pca_variance.csv` | Variance explained per PC |
| `fasting_pca_contributions.csv` | Variable contributions (squared loadings) |
| `fasting_pca_input_log_complete_cases.csv` | Log-transformed input matrix |
| `fasting_pca_missingness.csv` | Missing/non-positive row counts per compound |
| `fasting_pca_verification_checks.csv` | Data audit trail |

**PNG plots** (in `outputs/pca/plots/`):

| File | Description |
|------|-------------|
| `pca_scores_plot.png` | Score plot with ellipses and centroids |
| `pca_subject_trajectories.png` | Paired-subject movement in PCA space |
| `pca_scree_plot.png` | Variance explained summary |
| `pca_loadings_plot.png` | Loading directions coloured by bile acid class |
| `pca_score_distributions.png` | Box plots of PC scores by condition |
| `pca_top_loadings.png` | Top signed loadings for PC1–PC3 |
| `pca_variable_contributions.png` | Variable contribution percentages |
