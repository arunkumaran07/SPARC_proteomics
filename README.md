# sparcproteomics

Engine-agnostic downstream analysis of protein-level quantitative proteomics in R.
Point it at a DIA-NN, Spectronaut or MaxQuant protein table (or a generic
matrix, or an annotated Excel workbook) plus a sample sheet, choose your
options in one YAML file, and get differential-abundance tables, plots,
Reactome enrichment and, optionally, an AI-written interpretation.

> Status: early development (0.1.0). Private while testing.

## What you can choose

| Option | Choices |
|---|---|
| Species | human, mouse, rat, zebrafish, fly, worm, other (`list_species()`) |
| Input format | `diann`, `spectronaut`, `maxquant`, `generic`, `annotated_xlsx` |
| Control group | any level of your group column |
| Comparisons | every group vs control, all pairwise, or a list (`"DoseHigh vs DoseLow"`) |
| Number of groups | any (k ≥ 2); a global F-test is added for k > 2 |
| Filtering | valid values in any group / all groups / overall / none |
| Normalization | none, median, quantile, cyclic loess |
| Imputation | none, mixed (kNN + down-shift), kNN, down-shift (Perseus-style), minimum |
| Statistics | limma moderated t (recommended), Welch t, Student t |
| Covariates / blocking | any sample-sheet columns (limma) |
| AI interpretation | auto (only if an API key is set), on, off |

## Install

```r
install.packages("remotes")
remotes::install_github("arunkumaran07/SPARC_proteomics")   # needs access while private
# optional extras
BiocManager::install(c("impute", "clusterProfiler", "ReactomePA", "org.Hs.eg.db", "org.Mm.eg.db"))
install.packages(c("ggplot2", "readxl", "httr2"))
```

## Quick start

```r
library(sparcproteomics)

# 1. try it on simulated data
cfg <- simulate_dataset("demo", n_groups = 3)
res <- run_pipeline(cfg)
res$de$summary

# 2. your own data: answer questions in the console ...
create_config("my_study.yml")
run_pipeline("my_study.yml")

# ... or copy and edit the annotated template
file.copy(system.file("extdata", "example_config.yml", package = "sparcproteomics"), "my_study.yml")
```

The sample sheet is a CSV with one row per sample: a sample-ID column matching
the matrix columns (folders, file extensions and Spectronaut `[n]` prefixes are
ignored when matching), a group column, and any covariates.

## Design choices

- **Tests never see imputed values** by default. limma fits each protein on its
  observed values; imputation feeds only PCA and other methods that need a
  complete matrix.
- **limma is the default test.** On simulated data with 3 samples per group it
  recovered ~80% of true changes where unmoderated t-tests recovered almost
  none; the moderated statistic borrows variance information across proteins
  (Smyth 2004).
- **Per-comparison gating**: a protein is tested in a comparison only if both
  groups have at least `min_valid_per_group` values; gating happens before FDR
  correction.
- **Fail loudly**: unknown control group, missing design columns, repeated
  sample names and blank annotations stop the run or are logged and dropped,
  never silently mis-modelled.
- **AI is optional** and receives only summary statistics, never sample-level
  data. The key is read from `ANTHROPIC_API_KEY` and is never written out.

## Outputs (`output$dir`)

`DE_<comparison>.csv`, `DE_summary.csv`, `F_test.csv`, `Reactome_ORA_*.csv`,
`Reactome_GSEA_*.csv`, `PCA.png`, `Volcano_*.png`, `AI_interpretation.md`,
`config_used.yml`, `run_log.txt`, `session_info.txt`, `results.rds`
(the `SummarizedExperiment` and all results).

## Roadmap

QC report, tissue/contamination marker panels (muscle, plasma, CSF),
presence/absence testing, time-course trajectories, factorial interaction
tests, cross-validated classification, PDF report.

## Key references

- Smyth GK (2004) *Stat Appl Genet Mol Biol* 3:Article 3 — moderated t-statistics.
- Ritchie ME et al. (2015) *Nucleic Acids Res* 43:e47 — limma.
- Phipson B et al. (2016) *Ann Appl Stat* 10:946–963 — robust empirical Bayes.
- Tyanova S et al. (2016) *Nat Methods* 13:731–740 — Perseus (down-shift imputation).
- Yu G, He Q-Y (2016) *Mol BioSyst* 12:477–479 — ReactomePA.

## License

MIT
