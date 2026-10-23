#' Filter proteins by valid values
#'
#' @param se A `SummarizedExperiment` from [build_dataset()].
#' @param method `"any_group"` keeps proteins with at least `min_fraction`
#'   valid values in at least one group (keeps proteins present in only one
#'   condition); `"all_groups"` requires it in every group; `"overall"` across
#'   all samples; `"none"` removes only proteins with no values at all.
#' @param min_fraction Minimum fraction of valid (non-missing) values.
#' @param assay Assay to evaluate.
#' @return The filtered `SummarizedExperiment`. The pre-filter object is not
#'   kept; keep your own copy if you need it.
#' @export
pp_filter <- function(se, method = "any_group", min_fraction = 0.7, assay = "log2") {
  x <- SummarizedExperiment::assay(se, assay)
  valid <- !is.na(x)
  g <- se$group
  frac_by_group <- vapply(levels(g), function(l)
    rowMeans(valid[, g == l, drop = FALSE]), numeric(nrow(x)))
  if (is.null(dim(frac_by_group))) frac_by_group <- matrix(frac_by_group, nrow = 1)
  keep <- switch(method,
    any_group  = apply(frac_by_group >= min_fraction, 1, any),
    all_groups = apply(frac_by_group >= min_fraction, 1, all),
    overall    = rowMeans(valid) >= min_fraction,
    none       = rowSums(valid) > 0,
    stop("Unknown filter method '", method, "'.", call. = FALSE))
  log_msg("Filter (", method, if (method != "none") paste0(", >= ", min_fraction * 100, "% valid"),
          "): kept ", sum(keep), " of ", length(keep), " proteins.")
  out <- se[keep, ]
  S4Vectors::metadata(out)$filter <- list(method = method, min_fraction = min_fraction,
                                          n_before = length(keep), n_after = sum(keep))
  out
}

#' Normalize log2 intensities between samples
#'
#' @param se A `SummarizedExperiment` with a `log2` assay.
#' @param method `"none"` (data already normalized upstream, e.g. by the search
#'   engine), `"median"` (subtract each sample's median, add back the overall
#'   median), `"mode"` (see below), `"quantile"` or `"cyclicloess"` (both via
#'   [limma::normalizeBetweenArrays()]).
#'
#' Median normalization assumes most proteins do not change. When a large
#' fraction does (e.g. strong treatments, spike-in or mixed-species designs,
#' tissue-composition shifts), either give a `reference` set expected to be
#' stable (sample medians are then computed from those proteins only), or use
#' `"mode"`.
#'
#' `"mode"` compares samples in pairs: for each pair, the **most common**
#' per-protein log-ratio (the peak of its density) estimates how far apart the
#' two samples are, and the sample shifts are the least-squares fit to all
#' pairs. Proteins that change move away from that peak rather than dragging
#' it, so the method only assumes that unchanged proteins are the largest
#' single group, not a majority. The price is a noisier estimate: when changes
#' are balanced, median normalization is somewhat more precise.
#'
#' @param reference For `"median"`: optional protein IDs or gene symbols
#'   (matched against `protein_id` and `gene`) used to compute the
#'   normalization factors. All proteins are still normalized.
#' @return `se` with a new assay `norm`.
#' @export
pp_normalize <- function(se, method = "mode", reference = NULL) {
  x <- SummarizedExperiment::assay(se, "log2")
  ref_rows <- rep(TRUE, nrow(x))
  if (!is.null(reference) && method == "median") {
    rd <- SummarizedExperiment::rowData(se)
    ref_rows <- rd$protein_id %in% reference | (!is.na(rd$gene) & rd$gene %in% reference)
    if (sum(ref_rows) < 20)
      stop("Only ", sum(ref_rows), " reference proteins matched; need at least 20.", call. = FALSE)
    log_msg("Median normalization on ", sum(ref_rows), " reference proteins.")
  } else if (!is.null(reference)) {
    log_msg("NOTE: normalization reference is only used by method 'median'; ignored.")
  }
  y <- switch(method,
    none        = x,
    median      = sweep(x, 2, matrixStats::colMedians(x[ref_rows, , drop = FALSE], na.rm = TRUE), "-") +
                    stats::median(x[ref_rows, ], na.rm = TRUE),
    mode        = sweep(x, 2, .mode_shifts(x), "-"),
    quantile    = limma::normalizeBetweenArrays(x, method = "quantile"),
    cyclicloess = limma::normalizeBetweenArrays(x, method = "cyclicloess"),
    stop("Unknown normalization '", method, "'.", call. = FALSE))
  dimnames(y) <- dimnames(x)
  SummarizedExperiment::assay(se, "norm") <- y
  S4Vectors::metadata(se)$normalization <- method
  S4Vectors::metadata(se)$normalization_reference_n <- if (is.null(reference)) NA else sum(ref_rows)
  log_msg("Normalization: ", method, ".")
  se
}

#' Impute missing values
#'
#' Imputed values are intended for methods that cannot handle missing values
#' (PCA, clustering, heatmaps, classifiers). By default, differential testing
#' in [run_de()] uses the un-imputed `norm` assay so imputed values cannot
#' create significance.
#'
#' Methods:
#' \describe{
#'   \item{`none`}{No imputation; the `imputed` assay keeps the NAs.}
#'   \item{`downshift`}{Random draws from a normal distribution `shift` SDs
#'     below each **sample's** mean, with SD `width` x the sample SD
#'     (the Perseus approach; values assumed missing because they are below
#'     the detection limit).}
#'   \item{`mindet`}{Each sample's minimum observed value.}
#'   \item{`knn`}{k-nearest-neighbour imputation via [impute::impute.knn()]
#'     (values assumed missing at random).}
#'   \item{`mixed`}{Per group: a protein missing in at least
#'     `mnar_group_fraction` of a group's samples is treated as below the
#'     detection limit there and filled by `downshift`; remaining gaps are
#'     filled by kNN. Proteins still missing in >= 80% of samples after that
#'     get each sample's minimum.}
#' }
#'
#' @param se A `SummarizedExperiment` with a `norm` assay.
#' @param method See Details.
#' @param shift,width Down-shift parameters.
#' @param mnar_group_fraction Threshold for the `mixed` method.
#' @param k Number of neighbours for kNN.
#' @param seed Random seed.
#' @return `se` with assays `imputed` and `imputed_mask` (TRUE where a value
#'   was imputed).
#' @export
pp_impute <- function(se, method = "mixed", shift = 1.8, width = 0.3,
                      mnar_group_fraction = 0.6, k = 10, seed = 42) {
  x <- SummarizedExperiment::assay(se, "norm")
  miss <- is.na(x)
  set.seed(seed)
  downshift <- function(m, mask) {
    for (j in seq_len(ncol(m))) {
      idx <- which(mask[, j])
      if (!length(idx)) next
      mu <- mean(m[, j], na.rm = TRUE); s <- stats::sd(m[, j], na.rm = TRUE)
      m[idx, j] <- stats::rnorm(length(idx), mu - shift * s, width * s)
    }
    m
  }
  per_sample_min <- function(m, mask) {
    for (j in seq_len(ncol(m))) m[mask[, j], j] <- min(m[, j], na.rm = TRUE)
    m
  }
  knn <- function(m) {
    .need("impute", "for kNN imputation")
    rm <- rowMeans(is.na(m))
    if (any(rm > 0.5))
      log_msg("NOTE: ", sum(rm > 0.5), " proteins are >50% missing; impute.knn ",
              "fills these with the row mean, which understates their variance.")
    res <- utils::capture.output(out <- suppressWarnings(
      impute::impute.knn(m, k = min(k, ncol(m) - 1), rng.seed = seed)$data))
    out
  }

  y <- switch(method,
    none      = x,
    downshift = downshift(x, miss),
    mindet    = per_sample_min(x, miss),
    knn       = knn(x),
    mixed     = {
      g <- se$group
      mnar <- matrix(FALSE, nrow(x), ncol(x))
      for (l in levels(g)) {
        cols <- which(g == l)
        hit <- rowMeans(miss[, cols, drop = FALSE]) >= mnar_group_fraction
        mnar[hit, cols] <- miss[hit, cols]
      }
      y <- x
      ok <- rowMeans(miss & !mnar) < 0.8 & rowSums(!miss) >= 2
      if (any(ok) && any(miss[ok, ] & !mnar[ok, ])) {
        tmp <- x[ok, , drop = FALSE]
        tmp_k <- knn(tmp)
        fill <- miss[ok, , drop = FALSE] & !mnar[ok, , drop = FALSE]
        y_ok <- y[ok, , drop = FALSE]; y_ok[fill] <- tmp_k[fill]; y[ok, ] <- y_ok
      }
      y <- downshift(y, mnar)
      left <- is.na(y)
      if (any(left)) y <- per_sample_min(y, left)
      log_msg("Mixed imputation: ", sum(mnar), " values treated as below detection ",
              "(down-shift), ", sum(miss & !mnar & !left), " by kNN, ", sum(left),
              " by per-sample minimum.")
      y
    },
    stop("Unknown imputation method '", method, "'.", call. = FALSE))
  dimnames(y) <- dimnames(x)
  SummarizedExperiment::assay(se, "imputed") <- y
  SummarizedExperiment::assay(se, "imputed_mask") <- miss & !is.na(y)
  S4Vectors::metadata(se)$imputation <- list(method = method, shift = shift,
    width = width, mnar_group_fraction = mnar_group_fraction, k = k, seed = seed,
    n_missing = sum(miss), n_imputed = sum(miss & !is.na(y)))
  log_msg("Imputation (", method, "): ", sum(miss & !is.na(y)), " of ", sum(miss),
          " missing values filled.")
  se
}

# Location of the highest peak of a density estimate.
.density_mode <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 50) return(stats::median(v))
  d <- stats::density(v, n = 4096, bw = "nrd0")
  d$x[which.max(d$y)]
}

# Per-sample shifts for mode normalization. For pairs of samples, the mode of
# their per-protein log-ratios estimates the difference in their shifts; the
# shifts are then the least-squares fit to all pairwise differences (centred
# to sum to zero). Comparing samples pairwise avoids comparing a sample with a
# reference that contains it, which pulls the mode towards zero. Up to 61
# samples all pairs are used; above that each sample is paired with its next
# 30 neighbours (cyclically), which keeps every sample connected.
.mode_shifts <- function(x) {
  n <- ncol(x)
  use <- rowMeans(!is.na(x)) >= 0.5
  if (sum(use) < 50 || n < 2) {
    log_msg("NOTE: fewer than 50 well-observed proteins; mode normalization falls back to medians.")
    s <- matrixStats::colMedians(x, na.rm = TRUE)
    return(s - mean(s, na.rm = TRUE))
  }
  x <- x[use, , drop = FALSE]
  pr <- if (n <= 61) t(utils::combn(n, 2)) else
    do.call(rbind, lapply(seq_len(30), function(k) cbind(seq_len(n), (seq_len(n) + k - 1) %% n + 1)))
  r <- apply(pr, 1, function(ij) .density_mode(x[, ij[1]] - x[, ij[2]]))
  ok <- is.finite(r)
  D <- matrix(0, sum(ok), n)
  D[cbind(seq_len(sum(ok)), pr[ok, 1])] <- 1
  D[cbind(seq_len(sum(ok)), pr[ok, 2])] <- -1
  linked <- colSums(D != 0) > 0
  s <- rep(NA_real_, n)
  s[linked] <- qr.solve(rbind(D[, linked, drop = FALSE], 1), c(r[ok], 0))
  if (any(!linked)) {
    log_msg("NOTE: ", sum(!linked), " sample(s) share too few proteins with the others for ",
            "mode normalization; they are median-normalized instead.")
    med <- matrixStats::colMedians(x, na.rm = TRUE)
    s[!linked] <- med[!linked] - mean(med[linked], na.rm = TRUE) + mean(s[linked])
  }
  s
}

# Warns when the most common log2FC of a comparison is clearly not zero, the
# sign that a normalization assumption failed (many proteins changing in one
# direction). Returns the per-comparison modes invisibly.
.check_centering <- function(de, normalization, tolerance = 0.1) {
  modes <- vapply(de$tables, function(t) {
    v <- t$logFC[is.finite(t$logFC)]
    if (length(v)) .density_mode(v) else NA_real_
  }, numeric(1))
  off <- names(modes)[!is.na(modes) & abs(modes) > tolerance]
  if (length(off) && normalization != "mode")
    log_msg("NOTE: the most common log2FC is ", paste(sprintf("%.2f in %s", modes[off], off),
            collapse = ", "), ", not 0. With small groups this can be noise, but if many ",
            "proteins are expected to change in one direction, '", normalization,
            "' normalization shifts all ratios; consider preprocessing.normalization: mode ",
            "or a normalization_reference set (see benchmark/RESULTS.md).")
  invisible(modes)
}
