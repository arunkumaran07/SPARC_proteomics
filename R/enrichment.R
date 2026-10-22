#' Reactome pathway enrichment for one comparison
#'
#' Over-representation analysis of significant proteins against the background
#' of all proteins tested in that comparison, and gene-set enrichment analysis
#' on the moderated t-statistic. Gene symbols (the first symbol of each protein
#' group) are mapped to Entrez IDs with the species' `OrgDb` package.
#'
#' @param tab One table from `run_de()$tables`.
#' @param species A key from [list_species()].
#' @param min_gs_size,max_gs_size Gene-set size limits.
#' @param p_cutoff Adjusted p-value cutoff for reported pathways.
#' @return A list with `ora` and `gsea` data frames (possibly empty), and
#'   `status` describing what ran.
#' @export
run_enrichment <- function(tab, species, min_gs_size = 10, max_gs_size = 500,
                           p_cutoff = 0.05) {
  sp <- species_info(species)
  empty <- data.frame()
  if (is.na(sp$reactome))
    return(list(ora = empty, gsea = empty, status = "skipped: species not supported by Reactome"))
  for (p in c("clusterProfiler", "ReactomePA", sp$orgdb))
    if (!requireNamespace(p, quietly = TRUE))
      return(list(ora = empty, gsea = empty,
                  status = paste0("skipped: install ", p, " (BiocManager::install)")))
  orgdb <- getExportedValue(sp$orgdb, sp$orgdb)

  tab <- tab[!is.na(tab$P.Value) & !is.na(tab$gene), , drop = FALSE]
  tab$symbol <- sub("[;,].*$", "", tab$gene)
  map <- tryCatch(suppressMessages(suppressWarnings(
    clusterProfiler::bitr(unique(tab$symbol), fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = orgdb))), error = function(e) NULL)
  if (is.null(map) || nrow(map) < 20)
    return(list(ora = empty, gsea = empty, status = "skipped: too few symbols mapped to Entrez IDs"))
  map <- map[!duplicated(map$SYMBOL), ]
  tab <- merge(tab, map, by.x = "symbol", by.y = "SYMBOL")
  universe <- unique(tab$ENTREZID)
  sig <- unique(tab$ENTREZID[tab$significant])

  ora <- empty
  if (length(sig) >= 5) {
    r <- tryCatch(ReactomePA::enrichPathway(gene = sig, organism = sp$reactome,
      universe = universe, pvalueCutoff = p_cutoff, minGSSize = min_gs_size,
      maxGSSize = max_gs_size, readable = TRUE), error = function(e) NULL)
    if (!is.null(r)) ora <- as.data.frame(r)
  }
  ranked <- tab[order(-abs(tab$t)), ]
  ranked <- ranked[!duplicated(ranked$ENTREZID), ]
  gl <- sort(stats::setNames(ranked$t, ranked$ENTREZID), decreasing = TRUE)
  gsea <- empty
  r <- tryCatch(suppressWarnings(ReactomePA::gsePathway(gl, organism = sp$reactome,
    minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = p_cutoff,
    verbose = FALSE)), error = function(e) NULL)
  if (!is.null(r)) gsea <- as.data.frame(r)
  list(ora = ora, gsea = gsea,
       status = sprintf("ran: %d mapped, %d significant; %d ORA and %d GSEA pathways",
                        length(universe), length(sig), nrow(ora), nrow(gsea)))
}
