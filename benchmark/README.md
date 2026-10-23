# Benchmark

A check of `sparcproteomics` against real LC-MS data with known answers.
Results: [RESULTS.md](RESULTS.md).

## Dataset

- **Source:** PXD028735, from Van Puyvelde B *et al.* "A comprehensive LFQ benchmark
  dataset on modern day acquisition strategies in proteomics." *Sci Data* 9:126 (2022),
  [doi:10.1038/s41597-022-01216-6](https://doi.org/10.1038/s41597-022-01216-6).
- **Design:** commercial human (K562), yeast and *E. coli* digests mixed as
  A = 65 / 30 / 5 % and B = 65 / 15 / 20 % (w/w), 3 injections each.
  Expected log2(A/B): human 0, yeast +1, *E. coli* −2.
- **Subset used:** Q Exactive HF-X, data-dependent acquisition, quantified at precursor
  level with WOMBAT-P, as distributed in the
  [ProteoBench](https://github.com/Proteobench/ProteoBench) repository
  (Apache License 2.0; commit `2943465`).
- **Species labels:** taken from ProteoBench's MaxQuant result for the same runs, using
  species-unique precursors only. A protein group is kept when all its matched peptides
  agree on one species: 4,345 of 4,941 non-contaminant groups (3,023 human,
  1,076 yeast, 246 *E. coli*; 57,714 precursors).

`data/pxd028735_qehfx_dda_precursors.tsv.gz` is derived from the two ProteoBench files
named in `prepare_benchmark_data.R`, which regenerates it.

## Running it

From the repository root, with the package installed:

```r
Rscript benchmark/run_benchmark.R           # all 64 combinations -> RESULTS.md, results/
Rscript benchmark/run_benchmark.R --check   # also fail if results regress
```

It runs every combination of summarization (median polish, MaxLFQ, top3, sum),
normalization (median of all proteins, mode of ratios, median of a human reference set,
none) and test (limma with p + 1.5-fold cutoff, limma TREAT at 1.5-fold, Welch,
Student), 64 in all, starting from precursors, without imputation for testing. It also
scores the descriptive on/off (presence/absence) list. The
GitHub Actions workflow runs it with `--check` on every push.

## What it shows

- **limma vs t-tests.** With 3 vs 3 replicates, limma finds 78–88% of the true changes
  at ≤1.4% false discoveries. Welch finds 47–58% when data are normalized, and the
  unmoderated t-tests find nothing at all without normalization: their best p-value
  (~3e-4) cannot survive correction across ~4,300 proteins.
- **Summarization.** Median polish and MaxLFQ outperform top3 and sum on sensitivity
  and false discoveries in every normalization setting; MaxLFQ finds slightly more true
  changes with slightly more false ones.
- **Median normalization when many proteins change.** Here ~30% of proteins truly
  change, which violates the assumption that most do not. Normalizing on all proteins
  shifts every ratio by about −0.15 log2. The fold-change cutoff absorbs most of this
  (FDP 0.7%), but on adjusted p-values alone 17.7% of calls are false. Computing the
  factors from a stable reference set (`preprocessing.normalization_reference`)
  removes the shift (human bias 0.00) and raises sensitivity to 86%. Mode-of-ratios
  normalization (the package default) needs no reference list: 84% sensitivity,
  0.4% FDP, 1.9% FDP on adjusted p alone, human bias −0.04.
- **TREAT at 1.5-fold** produced no false discoveries in any normalization, but found
  10–45% of true changes (median polish), against 78–88% for the default rule.
- **Presence/absence** could not be evaluated: every protein here has values in both
  conditions.
- **Ratio compression.** The *E. coli* ratio comes out at about −1.1 to −1.5 instead of
  −2 in every setting, including no normalization, so it comes from the data (*E. coli*
  is 5% of sample A), not from the analysis. The CI thresholds allow for it.

## Limitations

One instrument and acquisition mode (Orbitrap DDA), one quantification tool, 3 vs 3
replicates, and a two-group design. It does not test DIA, timsTOF or SCIEX data, larger
cohorts, or covariates. A diaPASEF benchmark is the obvious next addition.
