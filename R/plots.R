#' PCA of samples
#'
#' Uses the imputed assay by default; with no imputation, only proteins
#' without missing values are used. Proteins are scaled to unit variance.
#'
#' @param se A processed `SummarizedExperiment`.
#' @param assay Assay to use.
#' @param color Sample-annotation column for colour.
#' @return A ggplot object.
#' @export
plot_pca <- function(se, assay = "imputed", color = "group") {
  .need("ggplot2", "for plotting")
  x <- SummarizedExperiment::assay(se, assay)
  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[matrixStats::rowVars(x) > 0, , drop = FALSE]
  if (nrow(x) < 3) stop("Too few complete proteins for PCA.", call. = FALSE)
  p <- stats::prcomp(t(x), scale. = TRUE)
  ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  df <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2],
                   colour = as.character(SummarizedExperiment::colData(se)[[color]]),
                   sample = colnames(se))
  ggplot2::ggplot(df, ggplot2::aes(.data$PC1, .data$PC2, colour = .data$colour)) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
                  colour = color,
                  title = sprintf("PCA (%d proteins, %s)", nrow(x), assay)) +
    ggplot2::theme_bw(base_size = 11)
}

#' Volcano plot for one comparison
#'
#' @param tab One table from `run_de()$tables`.
#' @param title Plot title.
#' @param adj_p_cutoff,log2fc_cutoff Threshold lines.
#' @param label_n Number of top significant proteins to label, split between up and down.
#' @param subtitle Optional subtitle; by default it states the cutoffs.
#' @return A ggplot object.
#' @export
plot_volcano <- function(tab, title = "", adj_p_cutoff = 0.05, log2fc_cutoff = 0.585,
                         label_n = 10, subtitle = NULL) {
  .need("ggplot2", "for plotting")
  d <- tab[!is.na(tab$P.Value), , drop = FALSE]
  d$class <- ifelse(d$significant, ifelse(d$logFC > 0, "up", "down"), "ns")
  d$label <- ifelse(is.na(d$gene), d$protein_id, sub("[;,].*$", "", d$gene))
  s <- d[d$significant, ]
  s <- s[order(s$P.Value), ]
  up <- s[s$logFC > 0, ]; dn <- s[s$logFC < 0, ]
  n_dn <- min(nrow(dn), label_n %/% 2)
  lab <- rbind(utils::head(up, label_n - n_dn), utils::head(dn, n_dn))
  # ggrepel (if installed) moves labels apart; otherwise overlapping labels are dropped
  labels <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = .data$label), size = 2.8,
                             max.overlaps = Inf, min.segment.length = 0,
                             segment.colour = "grey60", show.legend = FALSE)
  } else {
    ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$label), size = 2.8,
                       vjust = -0.6, check_overlap = TRUE, show.legend = FALSE)
  }
  ggplot2::ggplot(d, ggplot2::aes(.data$logFC, -log10(.data$P.Value), colour = .data$class)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * log2fc_cutoff, linetype = 2, colour = "grey50") +
    labels +
    ggplot2::scale_colour_manual(values = c(down = "#2166AC", ns = "grey70", up = "#B2182B")) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 p-value", colour = NULL,
                  title = title,
                  subtitle = subtitle %||% sprintf("%d significant (adj. p < %g, |log2FC| >= %g)",
                                                   sum(d$significant), adj_p_cutoff, log2fc_cutoff)) +
    ggplot2::theme_bw(base_size = 11)
}
