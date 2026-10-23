#' Summarize precursors or peptides to proteins
#'
#' Combines precursor- (or peptide-) level log2 intensities into one value per
#' protein and sample.
#'
#' \describe{
#'   \item{`median_polish`}{Tukey's median polish on the precursor x sample
#'     log2 matrix of each protein; the protein value is the overall effect plus
#'     the sample effect. Robust to outlying precursors and to precursors
#'     missing in some samples (the approach of MSstats' TMP summarization).}
#'   \item{`top3`}{For each protein, the three precursors with the highest
#'     mean log2 intensity across samples; the protein value is the log2 of the
#'     mean of their linear intensities in each sample (Silva et al. 2006).}
#'   \item{`sum`}{log2 of the summed linear intensities of all precursors.
#'     Simple, but a precursor missing in one sample lowers that sample's
#'     value, so it is sensitive to missing values.}
#' }
#' A protein gets `NA` in samples where none of its precursors was observed.
#'
#' @param x Numeric matrix of log2 intensities, precursors x samples.
#' @param protein Protein (group) identifier for each row of `x`.
#' @param method `"median_polish"`, `"top3"` or `"sum"`.
#' @param min_precursors Proteins with fewer precursors are dropped.
#' @return A list with `x` (proteins x samples log2 matrix) and
#'   `n_precursors` (named integer vector).
#' @export
#' @examples
#' x <- rbind(p1 = c(20, 21, 22), p2 = c(18, 19, NA), p3 = c(25, 25, 26))
#' summarize_precursors(x, c("A", "A", "B"))$x
summarize_precursors <- function(x, protein, method = c("median_polish", "top3", "sum"),
                                 min_precursors = 1) {
  method <- match.arg(method)
  stopifnot(nrow(x) == length(protein))
  protein <- as.character(protein)
  keep <- !is.na(protein) & nzchar(protein) & rowSums(!is.na(x)) > 0
  x <- x[keep, , drop = FALSE]; protein <- protein[keep]
  idx <- split(seq_len(nrow(x)), protein)
  n <- lengths(idx)
  idx <- idx[n >= min_precursors]
  n_unconverged <- 0
  one <- function(rows) {
    m <- x[rows, , drop = FALSE]
    seen <- colSums(!is.na(m)) > 0
    est <- if (nrow(m) == 1) m[1, ] else switch(method,
      median_polish = {
        r <- withCallingHandlers(
          stats::medpolish(m, na.rm = TRUE, trace.iter = FALSE, maxiter = 50),
          warning = function(w) if (grepl("converge", conditionMessage(w))) {
            n_unconverged <<- n_unconverged + 1
            invokeRestart("muffleWarning")
          })
        r$overall + r$col
      },
      top3 = {
        top <- utils::head(order(rowMeans(m, na.rm = TRUE), decreasing = TRUE), 3)
        log2(colMeans(2^m[top, , drop = FALSE], na.rm = TRUE))
      },
      sum = log2(colSums(2^m, na.rm = TRUE)))
    est[!seen] <- NA
    est
  }
  out <- t(vapply(idx, one, numeric(ncol(x))))
  if (ncol(x) == 1) out <- t(out)
  if (n_unconverged)
    log_msg("NOTE: median polish did not fully converge for ", n_unconverged,
            " protein(s) in 50 iterations; their estimates are still usable.")
  dimnames(out) <- list(names(idx), colnames(x))
  out[!is.finite(out)] <- NA
  list(x = out, n_precursors = n[names(idx)])
}

#' Read a precursor- or peptide-level table
#'
#' @param path File path.
#' @param format `"diann_report"`: the DIA-NN main report in long format
#'   (`report.tsv`, or `report.parquet` from DIA-NN 2.x, which needs the
#'   \pkg{arrow} or \pkg{nanoparquet} package), one row per precursor and run.
#'   `"peptide_generic"`: a wide table with one row per precursor or peptide,
#'   a protein column, and one numeric column per sample.
#' @param protein_column,precursor_column,gene_column Column names for
#'   `peptide_generic`. If `precursor_column` is `NULL`, rows are numbered.
#' @param quantity DIA-NN quantity column, default `"Precursor.Normalised"`.
#' @param qvalue DIA-NN filter: precursors are kept if `Q.Value` and (when
#'   present) `PG.Q.Value` are at or below this value.
#' @return A list with `intensities` (precursors x samples, linear or as given)
#'   and `features` (data frame with `precursor_id`, `protein_id`, `gene`).
#' @export
read_precursor_table <- function(path, format = c("diann_report", "peptide_generic"),
                                 protein_column = NULL, precursor_column = NULL,
                                 gene_column = NULL, quantity = "Precursor.Normalised",
                                 qvalue = 0.01) {
  format <- match.arg(format)
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (format == "diann_report") .read_diann_report(path, quantity, qvalue)
  else .read_peptide_generic(path, protein_column, precursor_column, gene_column)
}

.read_table_fast <- function(path, cols = NULL) {
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    df <- if (requireNamespace("arrow", quietly = TRUE)) {
      as.data.frame(arrow::read_parquet(path))
    } else if (requireNamespace("nanoparquet", quietly = TRUE)) {
      as.data.frame(nanoparquet::read_parquet(path))
    } else stop("Reading .parquet needs the 'arrow' or 'nanoparquet' package: ",
                "install.packages(\"nanoparquet\").", call. = FALSE)
    return(if (is.null(cols)) df else df[intersect(cols, names(df))])
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    hdr <- names(data.table::fread(path, nrows = 0))
    sel <- if (is.null(cols)) hdr else intersect(cols, hdr)
    return(as.data.frame(data.table::fread(path, select = sel, na.strings = c("", "NA", "NaN"),
                                           showProgress = FALSE)))
  }
  df <- utils::read.delim(path, check.names = FALSE, quote = "", stringsAsFactors = FALSE,
                          na.strings = c("", "NA", "NaN"))
  if (is.null(cols)) df else df[intersect(cols, names(df))]
}

.read_diann_report <- function(path, quantity, qvalue) {
  need <- c("Run", "Precursor.Id", "Protein.Group", "Genes", quantity, "Q.Value", "PG.Q.Value")
  df <- .read_table_fast(path, need)
  miss <- setdiff(c("Run", "Precursor.Id", "Protein.Group", quantity), names(df))
  if (length(miss))
    stop("Not a DIA-NN main report, or the quantity column is missing: ",
         paste(miss, collapse = ", "), call. = FALSE)
  n0 <- nrow(df)
  ok <- rep(TRUE, n0)
  if ("Q.Value" %in% names(df)) ok <- ok & !is.na(df$Q.Value) & df$Q.Value <= qvalue
  if ("PG.Q.Value" %in% names(df)) ok <- ok & !is.na(df$PG.Q.Value) & df$PG.Q.Value <= qvalue
  df <- df[ok, , drop = FALSE]
  log_msg("DIA-NN report: kept ", nrow(df), " of ", n0, " precursor rows at q <= ", qvalue,
          " (Q.Value", if ("PG.Q.Value" %in% names(df)) " and PG.Q.Value", ").")
  runs <- unique(df$Run)
  prec <- unique(df$Precursor.Id)
  m <- matrix(NA_real_, length(prec), length(runs), dimnames = list(prec, runs))
  m[cbind(match(df$Precursor.Id, prec), match(df$Run, runs))] <- as.numeric(df[[quantity]])
  first <- !duplicated(df$Precursor.Id)
  feats <- data.frame(precursor_id = df$Precursor.Id[first],
                      protein_id = df$Protein.Group[first],
                      gene = if ("Genes" %in% names(df)) df$Genes[first] else NA_character_,
                      stringsAsFactors = FALSE)
  feats <- feats[match(prec, feats$precursor_id), ]
  list(intensities = m, features = feats)
}

.read_peptide_generic <- function(path, protein_column, precursor_column, gene_column) {
  if (is.null(protein_column))
    stop("input.protein_column is required for format 'peptide_generic'.", call. = FALSE)
  df <- .read_table_fast(path)
  for (k in c(protein_column, precursor_column, gene_column))
    if (!k %in% names(df)) stop("Column '", k, "' not found.", call. = FALSE)
  id_cols <- c(protein_column, precursor_column, gene_column)
  num <- vapply(df, is.numeric, logical(1))
  samp <- setdiff(names(df)[num], id_cols)
  if (!length(samp)) stop("No numeric sample columns found.", call. = FALSE)
  m <- as.matrix(df[samp]); storage.mode(m) <- "double"
  pid <- if (is.null(precursor_column)) paste0("row", seq_len(nrow(df))) else df[[precursor_column]]
  rownames(m) <- make.unique(as.character(pid))
  feats <- data.frame(precursor_id = rownames(m), protein_id = as.character(df[[protein_column]]),
                      gene = if (is.null(gene_column)) NA_character_ else as.character(df[[gene_column]]),
                      stringsAsFactors = FALSE)
  log_msg("Peptide table: ", nrow(m), " rows x ", ncol(m), " samples.")
  list(intensities = m, features = feats)
}

# Reads a precursor-level input and summarizes it to proteins on the scale
# the rest of the pipeline expects (linear unless input.is_log2 is TRUE).
.read_and_summarize <- function(inp) {
  tab <- read_precursor_table(inp$matrix, inp$format, inp$protein_column,
                              inp$precursor_column, inp$gene_column,
                              inp$diann_quantity, inp$diann_qvalue)
  m <- tab$intensities
  m[is.nan(m)] <- NA
  if (isTRUE(inp$zero_is_missing)) m[!is.na(m) & m == 0] <- NA
  lg <- if (isTRUE(inp$is_log2)) m else log2(m)
  lg[!is.finite(lg)] <- NA
  s <- summarize_precursors(lg, tab$features$protein_id, inp$summarization, inp$min_precursors)
  genes <- tapply(tab$features$gene, tab$features$protein_id, function(g) g[!is.na(g)][1])
  log_msg("Summarized ", nrow(m), " precursors to ", nrow(s$x), " proteins (",
          inp$summarization, ", >= ", inp$min_precursors, " precursor(s)).")
  feats <- data.frame(protein_id = rownames(s$x), gene = unname(genes[rownames(s$x)]),
                      n_precursors = as.integer(s$n_precursors),
                      row.names = rownames(s$x), stringsAsFactors = FALSE)
  list(intensities = if (isTRUE(inp$is_log2)) s$x else 2^s$x, features = feats)
}

.precursor_formats <- c("diann_report", "peptide_generic")
