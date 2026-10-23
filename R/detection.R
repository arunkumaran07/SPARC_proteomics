#' Presence/absence (differential detection) testing
#'
#' Abundance tests cannot use a protein that is detected in one group and
#' absent from the other: it has no fold change and is gated out (or, if
#' imputed, gets an arbitrary one). This function tests the detection pattern
#' itself, on the matrix **before** valid-value filtering, using Fisher's exact
#' test on detected/not-detected counts per group (ported from the SPARC v1.44
#' pipeline).
#'
#' With small groups Fisher's test cannot reach significance: with 3 vs 3,
#' the most extreme pattern (3/3 vs 0/3) gives p = 0.1. The function therefore
#' also flags **on/off** proteins descriptively: detected in at least
#' `on_off_min_fraction` of one group's samples and in none of the other's.
#' These flags carry no p-value and should be read as candidates. Detection is
#' sensitive to per-sample depth, so check identification counts per sample
#' before interpreting them.
#'
#' @param se A `SummarizedExperiment` from [build_dataset()] (before filtering).
#' @param comparisons Data frame from [build_comparisons()].
#' @param adj_p_cutoff Adjusted p-value cutoff (BH within each comparison).
#' @param min_delta Minimum difference in detection rate for significance.
#' @param on_off_min_fraction Minimum detection rate in the "on" group for the
#'   descriptive on/off flag.
#' @param assay Assay whose missing values define "not detected".
#' @return A list with `tables` (one per comparison, sorted by p-value) and
#'   `summary`.
#' @export
run_detection <- function(se, comparisons, adj_p_cutoff = 0.05, min_delta = 0.5,
                          on_off_min_fraction = 2 / 3, assay = "log2") {
  x <- SummarizedExperiment::assay(se, assay)
  g <- se$group
  feats <- as.data.frame(SummarizedExperiment::rowData(se))
  tables <- lapply(seq_len(nrow(comparisons)), function(i) {
    a <- comparisons$numerator[i]; b <- comparisons$denominator[i]
    xa <- x[, g == a, drop = FALSE]; xb <- x[, g == b, drop = FALSE]
    na <- ncol(xa); nb <- ncol(xb)
    da <- rowSums(!is.na(xa)); db <- rowSums(!is.na(xb))
    seen <- da + db > 0
    p <- .fisher_counts(da, na, db, nb)
    p[!seen] <- NA
    fa <- da / na; fb <- db / nb
    tab <- data.frame(
      protein_id = feats$protein_id, gene = feats$gene,
      detected_a = da, n_a = na, detected_b = db, n_b = nb,
      fraction_a = fa, fraction_b = fb, delta_fraction = fa - fb,
      P.Value = p, adj.P.Val = stats::p.adjust(p, "BH"),
      mean_log2_a = suppressWarnings(rowMeans(xa, na.rm = TRUE)),
      mean_log2_b = suppressWarnings(rowMeans(xb, na.rm = TRUE)),
      row.names = rownames(x), stringsAsFactors = FALSE)
    tab$mean_log2_a[!is.finite(tab$mean_log2_a)] <- NA
    tab$mean_log2_b[!is.finite(tab$mean_log2_b)] <- NA
    tab$significant <- !is.na(tab$adj.P.Val) & tab$adj.P.Val < adj_p_cutoff &
      abs(tab$delta_fraction) >= min_delta
    tab$on_off <- ifelse(fa >= on_off_min_fraction & db == 0, paste0("only_in_", a),
                  ifelse(fb >= on_off_min_fraction & da == 0, paste0("only_in_", b), ""))
    # rename the per-group columns in one step, so a group name ending in
    # "_a" or "_b" cannot be renamed twice
    nm <- names(tab)
    is_a <- grepl("_a$", nm); is_b <- grepl("_b$", nm)
    nm[is_a] <- paste0(sub("_a$", "", nm[is_a]), "_", a)
    nm[is_b] <- paste0(sub("_b$", "", nm[is_b]), "_", b)
    names(tab) <- nm
    tab <- tab[seen, , drop = FALSE]
    tab[order(tab$P.Value, -abs(tab$delta_fraction)), ]
  })
  names(tables) <- comparisons$name

  min_p <- vapply(seq_len(nrow(comparisons)), function(i) {
    na <- sum(g == comparisons$numerator[i]); nb <- sum(g == comparisons$denominator[i])
    .fisher_counts(na, na, 0, nb)
  }, numeric(1))
  summary <- data.frame(
    comparison = comparisons$name,
    min_possible_p = min_p,
    significant = vapply(tables, function(t) sum(t$significant), integer(1)),
    on_off = vapply(tables, function(t) sum(nzchar(t$on_off)), integer(1)),
    row.names = NULL)
  for (i in which(min_p >= adj_p_cutoff))
    log_msg("NOTE: with these group sizes the detection test cannot reach p < ", adj_p_cutoff,
            " for ", comparisons$name[i], " (smallest possible p = ", signif(min_p[i], 2),
            "); on/off proteins are listed descriptively.")
  log_msg("Detection: ", paste(sprintf("%s %d significant, %d on/off", summary$comparison,
          summary$significant, summary$on_off), collapse = "; "))
  list(tables = tables, summary = summary,
       settings = list(adj_p_cutoff = adj_p_cutoff, min_delta = min_delta,
                       on_off_min_fraction = on_off_min_fraction))
}

# Two-sided Fisher exact p-values for detected/not-detected counts, computed
# once per distinct count pattern.
.fisher_counts <- function(da, na, db, nb) {
  key <- paste(da, db)
  uk <- unique(key)
  pv <- vapply(uk, function(k) {
    v <- as.integer(strsplit(k, " ")[[1]])
    stats::fisher.test(matrix(c(v[1], na - v[1], v[2], nb - v[2]), 2))$p.value
  }, numeric(1))
  unname(pv[key])
}
