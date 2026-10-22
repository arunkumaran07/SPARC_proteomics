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

# 1. try the bundled example (copied to ./sparc_example, results in ./sparc_example/results)
cfg <- example_data("diann", "sparc_example")
res <- run_pipeline(cfg)
res$de$summary

# 2. or simulate any design you like
cfg <- simulate_dataset("demo", n_groups = 4, n_per_group = 3)

# 3. your own data: answer questions in the console ...
create_config("my_study.yml")
run_pipeline("my_study.yml")

# ... or copy and edit the annotated template
file.copy(system.file("extdata", "example_config.yml", package = "sparcproteomics"), "my_study.yml")
```

The sample sheet is a CSV with one row per sample: a sample-ID column matching
the matrix columns (folders, file extensions and Spectronaut `[n]` prefixes are
ignored when matching), a group column, and any covariates.

## Example data

Two **simulated** datasets ship in `inst/extdata/` (no real samples). They use
real human gene symbols and UniProt accessions so annotation and Reactome
enrichment behave as on real data; `truth.csv` lists every planted change.

| `example_data()` | Contents |
|---|---|
| `"diann"` | DIA-NN `report.pg_matrix.tsv`, 1,800 proteins; Control, TreatmentA, TreatmentB x 5; sex and batch columns. TreatmentA raises electron-transport-chain proteins (GO:0022900), TreatmentB lowers cytoplasmic-translation proteins (GO:0002181), plus 40 random changes each. |
| `"annotated"` | The same proteins as a core-facility Excel workbook: annotations above the data, Control vs TreatmentA, each sample exported as linear and log2 under one name. |

On the DIA-NN example the default settings recover 96% (TreatmentA) and 93%
(TreatmentB) of the planted changes with 3.5% and 4.8% false discoveries.
`data-raw/make_example_data.R` regenerates both datasets.

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
