#' Build the list of comparisons
#'
#' @param levels Group levels; the first is the control.
#' @param spec `"all_vs_control"`, `"pairwise"`, or a character vector/list of
#'   comparisons written `"A - B"` or `"A vs B"` (A relative to B).
#' @param control Control level (defaults to `levels[1]`).
#' @return A data frame with `name`, `numerator` and `denominator`.
#' @export
#' @examples
#' build_comparisons(c("Control", "DoseLow", "DoseHigh"), "all_vs_control")
#' build_comparisons(c("Control", "DoseLow", "DoseHigh"), list("DoseHigh vs DoseLow"))
build_comparisons <- function(levels, spec = "all_vs_control", control = levels[1]) {
  spec <- unlist(spec)
  if (length(spec) == 1 && spec == "all_vs_control") {
    num <- setdiff(levels, control); den <- rep(control, length(num))
  } else if (length(spec) == 1 && spec == "pairwise") {
    pr <- utils::combn(levels, 2)
    num <- pr[2, ]; den <- pr[1, ]
  } else {
    parts <- strsplit(spec, "\\s+(vs\\.?|VS|-)\\s+")
    bad <- vapply(parts, length, integer(1)) != 2
    if (any(bad)) stop("Cannot parse comparison(s): ", paste(spec[bad], collapse = "; "),
                       ". Write them as 'A - B' or 'A vs B'.", call. = FALSE)
    num <- trimws(vapply(parts, `[`, "", 1)); den <- trimws(vapply(parts, `[`, "", 2))
    unk <- setdiff(c(num, den), levels)
    if (length(unk)) stop("Comparison names unknown group(s): ", paste(unk, collapse = ", "),
                          ". Groups are: ", paste(levels, collapse = ", "), call. = FALSE)
    if (any(num == den)) stop("A group cannot be compared with itself.", call. = FALSE)
  }
  data.frame(name = paste0(num, "_vs_", den), numerator = num, denominator = den,
             stringsAsFactors = FALSE)
}

#' Differential abundance testing
#'
#' Tests every configured comparison on the un-imputed `norm` assay (or the
#' `imputed` assay if `use_imputed = TRUE`). log2FC is numerator minus
#' denominator.
#'
#' \describe{
#'   \item{`limma`}{One linear model with a coefficient per group
#'     (`~ 0 + group`, plus any covariates), contrasts for each comparison,
#'     empirical-Bayes moderation (`trend`/`robust` configurable), and
#'     `duplicateCorrelation` when a blocking variable is given. With more than
#'     two groups, a moderated F-test across all groups is also reported.
#'     limma fits each protein on its observed values, so no imputation is
#'     needed.}
#'   \item{`welch`, `student`}{A two-sample t-test per comparison using only
#'     the samples of the two groups. Covariates and blocking are ignored.}
#' }
#' For each comparison, a protein's statistics are set to `NA` (before
#' multiple-testing correction) unless both groups have at least
#' `min_valid_per_group` observed values.
#'
#' Significance rule:
#' \describe{
#'   \item{`cutoff`}{adjusted p < `adj_p_cutoff` **and** |log2FC| >=
#'     `log2fc_cutoff` (the common double filter).}
#'   \item{`treat`}{limma only. TREAT (McCarthy & Smyth 2009) tests whether
#'     |log2FC| exceeds log2(`treat_fc`) and controls the FDR for that
#'     hypothesis; significant = adjusted p < `adj_p_cutoff`. It is stricter
#'     than a cutoff at the same fold change: a protein must be shown to exceed
#'     the threshold, not merely to have an estimate above it.}
#' }
#'
#' @param se A `SummarizedExperiment` after [pp_normalize()].
#' @param comparisons Data frame from [build_comparisons()].
#' @param method `"limma"`, `"welch"` or `"student"`.
#' @param covariates Sample-annotation columns to adjust for (limma only).
#' @param block Blocking column for repeated measures (limma only).
#' @param min_valid_per_group Minimum observed values per group.
#' @param p_adjust Multiple-testing method for [stats::p.adjust()], applied
#'   within each comparison.
#' @param adj_p_cutoff,log2fc_cutoff Significance thresholds.
#' @param trend,robust Passed to [limma::eBayes()].
#' @param global_f_test Report a moderated F-test across groups (limma, > 2
#'   groups).
#' @param use_imputed Test the imputed assay instead (not recommended).
#' @param significance `"cutoff"` or `"treat"`; see Details.
#' @param treat_fc Fold-change threshold for TREAT (linear scale, default 1.5).
#' @return A list: `tables` (one data frame per comparison), `summary`,
#'   `f_test` (or `NULL`), `method`, and `design` (limma only).
#' @export
run_de <- function(se, comparisons, method = "limma", covariates = NULL, block = NULL,
                   min_valid_per_group = 2, p_adjust = "BH", adj_p_cutoff = 0.05,
                   log2fc_cutoff = 0.585, trend = TRUE, robust = TRUE,
                   global_f_test = TRUE, use_imputed = FALSE,
                   significance = "cutoff", treat_fc = 1.5) {
  significance <- match.arg(significance, c("cutoff", "treat"))
  if (significance == "treat" && method != "limma") {
    log_msg("NOTE: TREAT is limma-only; using the fold-change cutoff for the ", method, " test.")
    significance <- "cutoff"
  }
  if (significance == "treat" && (!is.numeric(treat_fc) || treat_fc <= 1))
    stop("statistics.treat_fc must be a fold change greater than 1 (e.g. 1.5).", call. = FALSE)
  assay_name <- if (use_imputed) "imputed" else "norm"
  if (use_imputed) log_msg("WARNING: testing imputed values; imputed values can ",
                           "create apparent significance.")
  x <- SummarizedExperiment::assay(se, assay_name)
  g <- se$group
  feats <- as.data.frame(SummarizedExperiment::rowData(se))
  covariates <- unlist(covariates)

  res <- if (method == "limma") {
    .de_limma(x, g, comparisons, as.data.frame(SummarizedExperiment::colData(se)),
              covariates, block, trend, robust, global_f_test,
              treat_lfc = if (significance == "treat") log2(treat_fc) else NULL)
  } else {
    if (length(covariates) || !is.null(block))
      log_msg("NOTE: covariates/blocking are ignored by the ", method, " t-test; ",
              "use statistics.method = limma to adjust for them.")
    list(tables = lapply(seq_len(nrow(comparisons)), function(i)
      .de_ttest(x, g, comparisons$numerator[i], comparisons$denominator[i],
                equal_var = method == "student")), f_test = NULL, design = NULL)
  }
  names(res$tables) <- comparisons$name

  for (i in seq_len(nrow(comparisons))) {
    a <- comparisons$numerator[i]; b <- comparisons$denominator[i]
    tt <- res$tables[[i]]
    na <- rowSums(!is.na(x[, g == a, drop = FALSE]))
    nb <- rowSums(!is.na(x[, g == b, drop = FALSE]))
    gate <- na < min_valid_per_group | nb < min_valid_per_group
    tt[gate, c("logFC", "t", "P.Value")] <- NA
    tt$adj.P.Val <- stats::p.adjust(tt$P.Value, method = p_adjust)
    tt$significant <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < adj_p_cutoff &
      (significance == "treat" | abs(tt$logFC) >= log2fc_cutoff)
    tt$direction <- ifelse(is.na(tt$logFC), NA, ifelse(tt$logFC > 0, "up", "down"))
    tt[[paste0("n_valid_", a)]] <- na
    tt[[paste0("n_valid_", b)]] <- nb
    tt[[paste0("mean_", a)]] <- rowMeans(x[, g == a, drop = FALSE], na.rm = TRUE)
    tt[[paste0("mean_", b)]] <- rowMeans(x[, g == b, drop = FALSE], na.rm = TRUE)
    tt <- cbind(feats[rownames(tt), c("protein_id", "gene"), drop = FALSE], tt)
    res$tables[[i]] <- tt[order(tt$P.Value, na.last = TRUE), ]
  }

  res$summary <- data.frame(
    comparison = comparisons$name,
    numerator = comparisons$numerator, denominator = comparisons$denominator,
    n_numerator = as.integer(table(g)[comparisons$numerator]),
    n_denominator = as.integer(table(g)[comparisons$denominator]),
    tested = vapply(res$tables, function(t) sum(!is.na(t$P.Value)), integer(1)),
    significant = vapply(res$tables, function(t) sum(t$significant), integer(1)),
    up = vapply(res$tables, function(t) sum(t$significant & t$logFC > 0, na.rm = TRUE), integer(1)),
    down = vapply(res$tables, function(t) sum(t$significant & t$logFC < 0, na.rm = TRUE), integer(1)),
    min_adj_p = vapply(res$tables, function(t)
      if (all(is.na(t$adj.P.Val))) NA_real_ else min(t$adj.P.Val, na.rm = TRUE), numeric(1)),
    row.names = NULL)
  res$method <- method
  res$settings <- list(p_adjust = p_adjust, adj_p_cutoff = adj_p_cutoff,
                       log2fc_cutoff = log2fc_cutoff, significance = significance,
                       treat_fc = if (significance == "treat") treat_fc else NULL,
                       min_valid_per_group = min_valid_per_group, assay = assay_name,
                       covariates = covariates, block = block)
  log_msg("DE (", method, "): ", paste(sprintf("%s %d sig (%d up, %d down)",
          res$summary$comparison, res$summary$significant, res$summary$up,
          res$summary$down), collapse = "; "))
  res
}

.de_limma <- function(x, g, comparisons, cd, covariates, block, trend, robust, global_f_test,
                      treat_lfc = NULL) {
  safe <- paste0("G", seq_along(levels(g)))
  names(safe) <- levels(g)
  design <- stats::model.matrix(~ 0 + g)
  colnames(design) <- safe
  for (cv in covariates) {
    mm <- stats::model.matrix(stats::as.formula(paste("~", "cov")),
                              data = data.frame(cov = cd[[cv]]))[, -1, drop = FALSE]
    colnames(mm) <- make.names(paste0(cv, "_", sub("^cov", "", colnames(mm))))
    design <- cbind(design, mm)
  }
  if (qr(design)$rank < ncol(design))
    stop("The design is not estimable (a covariate is confounded with the groups). ",
         "Remove it from design.covariates.", call. = FALSE)
  rdf <- nrow(design) - ncol(design)
  log_msg("limma design: ", ncol(design), " coefficients, residual df = ", rdf, ".")
  if (rdf < 2) stop("Too few residual degrees of freedom (", rdf, ").", call. = FALSE)

  cor <- NULL
  if (!is.null(block)) {
    b <- cd[[block]]
    if (length(unique(b)) < length(b)) {
      dc <- tryCatch(limma::duplicateCorrelation(x, design, block = b),
                     error = function(e) NULL)
      if (!is.null(dc) && is.finite(dc$consensus) && dc$consensus > 0) {
        cor <- dc$consensus
        log_msg(sprintf("duplicateCorrelation on '%s': consensus = %.3f", block, cor))
      } else log_msg("Blocking on '", block, "' dropped (non-positive consensus correlation).")
    } else log_msg("Blocking variable '", block, "' has one sample per level; ignored.")
  }
  # "Partial NA coefficients" is expected when a protein is absent from a whole
  # group; those proteins are gated out per comparison in run_de().
  fit <- withCallingHandlers(
    if (is.null(cor)) limma::lmFit(x, design) else
      limma::lmFit(x, design, block = cd[[block]], correlation = cor),
    warning = function(w) if (grepl("Partial NA coefficients", conditionMessage(w))) {
      log_msg("NOTE: ", trimws(conditionMessage(w)), " (protein absent from a group; ",
              "affected comparisons report NA for it).")
      invokeRestart("muffleWarning")
    })

  cm <- matrix(0, ncol(design), nrow(comparisons),
               dimnames = list(colnames(design), comparisons$name))
  for (i in seq_len(nrow(comparisons))) {
    cm[safe[[comparisons$numerator[i]]], i] <- 1
    cm[safe[[comparisons$denominator[i]]], i] <- -1
  }
  cfit <- limma::contrasts.fit(fit, cm)
  if (is.null(treat_lfc)) {
    fit2 <- limma::eBayes(cfit, trend = trend, robust = robust)
    tables <- lapply(colnames(cm), function(ct) {
      tt <- limma::topTable(fit2, coef = ct, number = Inf, sort.by = "none")
      tt[rownames(x), c("logFC", "AveExpr", "t", "P.Value")]
    })
  } else {
    # TREAT: tests H0 |log2FC| <= treat_lfc, so p-values already account for
    # the fold-change threshold (McCarthy & Smyth 2009)
    fit2 <- limma::treat(cfit, lfc = treat_lfc, trend = trend, robust = robust)
    tables <- lapply(colnames(cm), function(ct) {
      tt <- limma::topTreat(fit2, coef = ct, number = Inf, sort.by = "none")
      tt[rownames(x), c("logFC", "AveExpr", "t", "P.Value")]
    })
  }

  f_test <- NULL
  if (global_f_test && nlevels(g) > 2) {
    cmf <- matrix(0, ncol(design), nlevels(g) - 1,
                  dimnames = list(colnames(design), paste0(levels(g)[-1], "_vs_", levels(g)[1])))
    for (j in seq_len(ncol(cmf))) { cmf[safe[j + 1], j] <- 1; cmf[safe[1], j] <- -1 }
    ff <- limma::eBayes(limma::contrasts.fit(fit, cmf), trend = trend, robust = robust)
    f_test <- data.frame(protein = rownames(x), F = ff$F, P.Value = ff$F.p.value,
                         adj.P.Val = stats::p.adjust(ff$F.p.value, "BH"), row.names = NULL)
    f_test <- f_test[order(f_test$P.Value, na.last = TRUE), ]
    log_msg("Global F-test across ", nlevels(g), " groups: ",
            sum(f_test$adj.P.Val < 0.05, na.rm = TRUE), " proteins at BH < 0.05.")
  }
  list(tables = tables, f_test = f_test, design = design)
}

.de_ttest <- function(x, g, a, b, equal_var = FALSE) {
  xa <- x[, g == a, drop = FALSE]; xb <- x[, g == b, drop = FALSE]
  na <- rowSums(!is.na(xa)); nb <- rowSums(!is.na(xb))
  ma <- rowMeans(xa, na.rm = TRUE); mb <- rowMeans(xb, na.rm = TRUE)
  va <- matrixStats::rowVars(xa, na.rm = TRUE); vb <- matrixStats::rowVars(xb, na.rm = TRUE)
  if (equal_var) {
    df <- na + nb - 2
    sp2 <- ((na - 1) * va + (nb - 1) * vb) / df
    se <- sqrt(sp2 * (1 / na + 1 / nb))
  } else {
    se2 <- va / na + vb / nb
    se <- sqrt(se2)
    df <- se2^2 / ((va / na)^2 / pmax(na - 1, 1) + (vb / nb)^2 / pmax(nb - 1, 1))
  }
  tstat <- (ma - mb) / se
  p <- 2 * stats::pt(-abs(tstat), df)
  bad <- na < 2 | nb < 2 | !is.finite(se) | se <= 0
  tstat[bad] <- NA; p[bad] <- NA
  data.frame(logFC = ma - mb, AveExpr = rowMeans(cbind(xa, xb), na.rm = TRUE),
             t = tstat, P.Value = p, row.names = rownames(x))
}

#' Describe the significance rule in words
#'
#' @param settings The `settings` element of a [run_de()] result.
#' @return A short character string, e.g. for figure subtitles and reports.
#' @export
significance_label <- function(settings) {
  if (identical(settings$significance, "treat"))
    sprintf("TREAT, |FC| > %g, %s-adjusted p < %g", settings$treat_fc, settings$p_adjust,
            settings$adj_p_cutoff)
  else sprintf("%s-adjusted p < %g and |log2FC| >= %g", settings$p_adjust %||% "BH",
               settings$adj_p_cutoff, settings$log2fc_cutoff)
}
