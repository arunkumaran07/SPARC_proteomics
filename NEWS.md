# sparcproteomics 0.2.0

## New options

* **Mode-of-ratios normalization is the new default** (`preprocessing.normalization:
  mode`). Samples are compared in pairs; the most common per-protein log-ratio of each
  pair estimates how far apart they are, and the sample shifts are fitted to all pairs.
  Proteins that change in one direction no longer bias the result. On the PXD028735
  benchmark it found 84% of true changes (median normalization: 78%) with 0.4% false
  discoveries (0.7%), and cut false calls on adjusted p alone from 17.7% to 1.9%. In
  simulations with balanced changes it found as many true changes as median with fewer
  false ones. Median normalization remains available (`normalization: median`).
* **MaxLFQ summarization** (`input.summarization: maxlfq`) for precursor-level input:
  an independent implementation of the published algorithm (Cox et al. 2014). Benchmark
  results are close to median polish, which stays the default.
* **TREAT** (`statistics.significance: treat`, `statistics.treat_fc: 1.5`): limma's
  test of whether a change exceeds the fold-change threshold, with FDR control for that
  hypothesis. Much stricter than the default p + 1.5-fold cutoff (benchmark: 0 false
  discoveries, 35% of true changes found with mode normalization).
* **Presence/absence testing** (`detection` section, `run_detection()`), ported from the
  SPARC v1.44 pipeline: Fisher's exact test on detection counts, run before valid-value
  filtering, plus a descriptive on/off list for small groups where Fisher's test cannot
  reach significance. Results in `Detection_<comparison>.csv`.

## Other changes

* A note in the run log when, under median or other normalization, the most common
  log2FC of a comparison is away from zero.
* License changed to Apache-2.0; NOTICE and CITATION.cff added.
* `significance_label()` describes the significance rule; volcano subtitles and the AI
  prompt use it.
* The benchmark now covers 64 combinations (4 summarizations x 4 normalizations x
  4 tests/rules) and scores on/off calls.

# sparcproteomics 0.1.0

* First version: configurable species, groups, control, comparisons, imputation and
  statistics; DIA-NN, Spectronaut, MaxQuant, generic and annotated-workbook input;
  DIA-NN precursor reports and peptide tables with median polish, top3 or sum
  summarization; reference-set median normalization; optional Reactome enrichment and
  AI interpretation; simulated example data; real-data benchmark on PXD028735.
