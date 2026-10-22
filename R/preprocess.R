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
#'   median), `"quantile"` or `"cyclicloess"` (both via
#'   [limma::normalizeBetweenArrays()]).
#' @return `se` with a new assay `norm`.
#' @export
pp_normalize <- function(se, method = "median") {
  x <- SummarizedExperiment::assay(se, "log2")
  y <- switch(method,
    none        = x,
    median      = sweep(x, 2, matrixStats::colMedians(x, na.rm = TRUE), "-") +
                    stats::median(x, na.rm = TRUE),
    quantile    = limma::normalizeBetweenArrays(x, method = "quantile"),
    cyclicloess = limma::normalizeBetweenArrays(x, method = "cyclicloess"),
    stop("Unknown normalization '", method, "'.", call. = FALSE))
  dimnames(y) <- dimnames(x)
  SummarizedExperiment::assay(se, "norm") <- y
  S4Vectors::metadata(se)$normalization <- method
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
