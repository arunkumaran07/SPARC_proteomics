# Benchmark results

Generated 2026-10-23 with sparcproteomics 0.2.0. Rerun with `Rscript benchmark/run_benchmark.R`.

**Dataset:** PXD028735 (Van Puyvelde et al. 2022), Q Exactive HF-X DDA, condition A vs B,
3 replicates each; precursor-level quantification from ProteoBench (WOMBAT-P).
1076 yeast, 245 *E. coli* and 3023 human proteins tested with default settings.
Expected log2(A/B): human 0, yeast +1, *E. coli* -2.

**Defaults scored:** median_polish summarization, mode of ratios normalization, limma, p + 1.5-fold cutoff.

**Scoring:** significant as reported by the pipeline: BH-adjusted p < 0.05 and |log2FC| >= 0.585
(1.5-fold), or for TREAT, BH-adjusted p < 0.05 for the hypothesis |FC| > 1.5. A true positive is
a yeast or *E. coli* protein called significant in the expected direction; any other
significant call is false. Sensitivity = true positives / yeast + *E. coli* proteins tested.
FDP = false / all significant calls. "FDP, adj. p only" ignores the fold-change requirement.
Bias = median log2FC minus the expected value. No imputation is used for testing.

## Significance rule and test (median_polish summarization)

| Normalization | Test and rule | Significant | Sensitivity | FDP | FDP, adj. p only |
|---|---|---|---|---|---|
| median (all proteins) | limma, p + 1.5-fold cutoff | 1037 | 78.0% | 0.7% | 17.7% |
| median (all proteins) | limma TREAT, 1.5-fold | 279 | 21.1% | 0.0% | 0.0% |
| median (all proteins) | Welch, p + 1.5-fold cutoff | 618 | 46.8% | 0.0% | 12.4% |
| median (all proteins) | Student, p + 1.5-fold cutoff | 913 | 69.1% | 0.0% | 17.2% |
| mode of ratios | limma, p + 1.5-fold cutoff | 1120 | 84.4% | 0.4% | 1.9% |
| mode of ratios | limma TREAT, 1.5-fold | 463 | 35.0% | 0.0% | 0.0% |
| mode of ratios | Welch, p + 1.5-fold cutoff | 699 | 52.9% | 0.0% | 1.0% |
| mode of ratios | Student, p + 1.5-fold cutoff | 958 | 72.5% | 0.0% | 1.2% |

## Summarization and normalization (limma, p + 1.5-fold cutoff)

| Summarization | Normalization | Sensitivity | FDP | FDP, adj. p only | Bias human | Bias yeast | Bias E. coli |
|---|---|---|---|---|---|---|---|
| median_polish | median (all proteins) | 78.0% | 0.7% | 17.7% | -0.15 | -0.27 | 0.49 |
| median_polish | mode of ratios | 84.4% | 0.4% | 1.9% | -0.04 | -0.17 | 0.59 |
| median_polish | median (human reference) | 86.1% | 0.4% | 2.6% | 0.00 | -0.13 | 0.64 |
| median_polish | none | 88.4% | 1.4% | 1.3% | 0.06 | -0.07 | 0.69 |
| maxlfq | median (all proteins) | 82.3% | 1.3% | 16.3% | -0.12 | -0.24 | 0.48 |
| maxlfq | mode of ratios | 87.4% | 0.6% | 2.2% | -0.02 | -0.14 | 0.58 |
| maxlfq | median (human reference) | 87.7% | 0.7% | 2.7% | -0.01 | -0.12 | 0.59 |
| maxlfq | none | 89.6% | 1.3% | 1.2% | 0.07 | -0.05 | 0.66 |
| top3 | median (all proteins) | 52.5% | 9.5% | 12.2% | -0.13 | -0.28 | 0.65 |
| top3 | mode of ratios | 57.7% | 11.0% | 14.3% | -0.01 | -0.16 | 0.77 |
| top3 | median (human reference) | 56.5% | 10.4% | 13.6% | -0.04 | -0.19 | 0.75 |
| top3 | none | 54.3% | 9.6% | 9.6% | 0.09 | -0.06 | 0.87 |
| sum | median (all proteins) | 63.7% | 5.5% | 14.6% | -0.13 | -0.22 | 0.45 |
| sum | mode of ratios | 70.3% | 5.5% | 10.6% | -0.00 | -0.10 | 0.58 |
| sum | median (human reference) | 69.6% | 5.4% | 10.1% | -0.02 | -0.12 | 0.56 |
| sum | none | 71.5% | 5.3% | 5.3% | 0.07 | -0.02 | 0.66 |

## Presence/absence (on/off) proteins

With 3 vs 3 samples Fisher's exact test cannot reach significance (smallest possible p = 0.1),
so on/off proteins are listed descriptively. Correct = yeast present only in A, or *E. coli*
present only in B.

| Rule | Listed | Correct | Human (false) | Wrong direction |
|---|---|---|---|---|
| >= 2 of 3 in one group, 0 in the other | 0 | 0 | 0 | 0 |
| 3 of 3 in one group, 0 in the other | 0 | 0 | 0 | 0 |

**Not evaluable on this dataset:** every protein is quantified in at least one run of each condition, so there are no on/off proteins to find. The module is covered by unit tests; benchmarking it needs a dataset with proteins absent from a whole condition.

All 64 combinations, with precision (MAD) per species, are in `results/metrics.csv`.

![Default settings](results/default_settings.png)
