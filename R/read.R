#' Read a sample sheet
#'
#' A sample sheet has one row per sample and one column per annotation: a
#' sample-ID column matching the protein table's sample columns, a group
#' column, and any covariates or blocking variables.
#'
#' @param path CSV, TSV/TXT or XLSX file.
#' @return A data frame with all columns as character.
#' @export
read_sample_sheet <- function(path) {
  if (!file.exists(path)) stop("Sample sheet not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  df <- if (ext %in% c("xlsx", "xls")) {
    .need("readxl", "to read Excel sample sheets")
    as.data.frame(readxl::read_excel(path, col_types = "text"))
  } else {
    utils::read.delim(path, sep = if (ext == "csv") "," else "\t",
                      colClasses = "character", check.names = FALSE,
                      na.strings = c("", "NA"))
  }
  names(df) <- trimws(names(df))
  df[] <- lapply(df, function(x) trimws(as.character(x)))
  df
}

#' Read a protein-group table
#'
#' @param path Path to the table.
#' @param format One of `"diann"` (`report.pg_matrix.tsv`), `"spectronaut"`
#'   (protein pivot report), `"maxquant"` (`proteinGroups.txt`), `"generic"`
#'   (one ID column plus numeric sample columns) or `"annotated_xlsx"` (see
#'   [read_annotated_xlsx()]).
#' @param intensity MaxQuant only: `"auto"` (LFQ if present, else raw
#'   intensity), `"lfq"`, `"intensity"` or `"ibaq"`.
#' @param id_column,gene_column Generic format only: names of the protein-ID
#'   and gene-symbol columns. If `id_column` is `NULL`, the first non-numeric
#'   column is used.
#' @param annotation_rows Passed to [read_annotated_xlsx()].
#' @return A list with `intensities` (proteins x samples numeric matrix),
#'   `features` (data frame with `protein_id`, `gene` and any description), and
#'   for `annotated_xlsx` also `samples`.
#' @export
read_protein_table <- function(path, format = c("diann", "spectronaut", "maxquant",
                                                "generic", "annotated_xlsx"),
                               intensity = "auto", id_column = NULL,
                               gene_column = NULL, annotation_rows = NULL) {
  format <- match.arg(format)
  if (!file.exists(path)) stop("Protein table not found: ", path, call. = FALSE)
  if (format == "annotated_xlsx")
    return(read_annotated_xlsx(path, annotation_rows = annotation_rows))

  tab <- utils::read.delim(path, check.names = FALSE, quote = "",
                           stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"),
                           comment.char = "")
  out <- switch(format,
    diann       = .parse_diann(tab),
    spectronaut = .parse_spectronaut(tab),
    maxquant    = .parse_maxquant(tab, intensity),
    generic     = .parse_generic(tab, id_column, gene_column))
  log_msg("Read ", nrow(out$intensities), " protein groups x ",
          ncol(out$intensities), " samples (", format, ").")
  out
}

.to_numeric_matrix <- function(df) {
  m <- suppressWarnings(as.matrix(data.frame(lapply(df, as.numeric),
                                             check.names = FALSE)))
  m[is.nan(m)] <- NA
  m
}

.finish <- function(m, ids, genes, desc = NULL) {
  ids <- as.character(ids)
  ids[is.na(ids) | !nzchar(ids)] <- paste0("unnamed_", which(is.na(ids) | !nzchar(ids)))
  genes <- as.character(genes %||% rep(NA_character_, length(ids)))
  genes[!nzchar(genes)] <- NA
  rownames(m) <- make.unique(ids)
  feats <- data.frame(protein_id = ids, gene = genes,
                      row.names = rownames(m), stringsAsFactors = FALSE)
  if (!is.null(desc)) feats$description <- as.character(desc)
  list(intensities = m, features = feats)
}

.parse_diann <- function(tab) {
  anno <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
            "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")
  if (!"Protein.Group" %in% names(tab))
    stop("This does not look like a DIA-NN pg_matrix (no 'Protein.Group' column).",
         call. = FALSE)
  samp <- setdiff(names(tab), anno)
  .finish(.to_numeric_matrix(tab[samp]), tab$Protein.Group, tab$Genes,
          tab$First.Protein.Description)
}

.parse_spectronaut <- function(tab) {
  if (!"PG.ProteinGroups" %in% names(tab))
    stop("This does not look like a Spectronaut protein pivot report ",
         "(no 'PG.ProteinGroups' column).", call. = FALSE)
  samp <- grep("PG\\.Quantity$", names(tab), value = TRUE)
  if (!length(samp)) samp <- setdiff(names(tab), grep("^PG\\.", names(tab), value = TRUE))
  m <- .to_numeric_matrix(tab[samp])
  colnames(m) <- sub("\\.PG\\.Quantity$", "", sub("^\\[[0-9]+\\]\\s*", "", samp))
  .finish(m, tab$PG.ProteinGroups, tab$PG.Genes, tab$PG.ProteinDescriptions)
}

.parse_maxquant <- function(tab, intensity) {
  if (!"Protein IDs" %in% names(tab))
    stop("This does not look like a MaxQuant proteinGroups.txt (no 'Protein IDs').",
         call. = FALSE)
  for (flag in c("Reverse", "Potential contaminant", "Only identified by site")) {
    if (flag %in% names(tab)) {
      hit <- !is.na(tab[[flag]]) & tab[[flag]] == "+"
      if (any(hit)) log_msg("MaxQuant: removed ", sum(hit), " '", flag, "' rows.")
      tab <- tab[!hit, , drop = FALSE]
    }
  }
  prefixes <- c(lfq = "^LFQ intensity ", intensity = "^Intensity ", ibaq = "^iBAQ ")
  if (identical(intensity, "auto"))
    intensity <- if (any(grepl(prefixes["lfq"], names(tab)))) "lfq" else "intensity"
  if (!intensity %in% names(prefixes))
    stop("MaxQuant intensity must be auto, lfq, intensity or ibaq.", call. = FALSE)
  samp <- grep(prefixes[[intensity]], names(tab), value = TRUE)
  if (!length(samp)) stop("No '", sub("\\^", "", prefixes[[intensity]]),
                          "<sample>' columns found.", call. = FALSE)
  log_msg("MaxQuant: using ", intensity, " columns.")
  m <- .to_numeric_matrix(tab[samp])
  colnames(m) <- sub(prefixes[[intensity]], "", samp)
  .finish(m, tab[["Protein IDs"]], tab[["Gene names"]], tab[["Protein names"]])
}

.parse_generic <- function(tab, id_column, gene_column) {
  num <- vapply(tab, function(x) {
    v <- suppressWarnings(as.numeric(x))
    mean(!is.na(v) | is.na(x)) > 0.9
  }, logical(1))
  if (is.null(id_column)) {
    if (!any(!num)) stop("No non-numeric ID column found; set input.id_column.",
                         call. = FALSE)
    id_column <- names(tab)[!num][1]
  }
  for (k in c(id_column, gene_column))
    if (!k %in% names(tab)) stop("Column '", k, "' not found.", call. = FALSE)
  samp <- setdiff(names(tab)[num], c(id_column, gene_column))
  .finish(.to_numeric_matrix(tab[samp]), tab[[id_column]],
          if (!is.null(gene_column)) tab[[gene_column]] else tab[[id_column]])
}

#' Read an Excel workbook with sample annotations above the data
#'
#' The layout used by many core-facility exports: the first row holds sample
#' names, the next rows hold annotations (group, sex, timepoint ...), with
#' labels in the first column, and the remaining rows hold one protein or
#' gene per row.
#'
#' @param path XLSX file.
#' @param annotation_rows Labels (first-column values) of the annotation rows.
#'   If `NULL`, they are detected as the contiguous rows under the header that
#'   are mostly non-numeric or have few distinct values. Giving them
#'   explicitly is safer.
#' @param sheet Sheet name or number.
#' @return A list with `intensities`, `features` and `samples` (data frame
#'   with `sample_id` plus one column per annotation row).
#' @export
read_annotated_xlsx <- function(path, annotation_rows = NULL, sheet = 1) {
  .need("readxl", "to read annotated Excel workbooks")
  raw <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_names = FALSE, .name_repair = "minimal",
                                          col_types = "text"))
  labels <- trimws(as.character(raw[[1]]))
  ids <- trimws(as.character(unlist(raw[1, -1])))
  if (any(is.na(ids) | !nzchar(ids)))
    stop("Blank sample name(s) in header row at column(s): ",
         paste(which(is.na(ids) | !nzchar(ids)) + 1, collapse = ", "), call. = FALSE)

  if (is.null(annotation_rows)) {
    i <- 2
    while (i <= nrow(raw)) {
      v <- as.character(unlist(raw[i, -1]))
      present <- !is.na(v) & nzchar(v)
      if (mean(present) < 0.5) break
      num <- suppressWarnings(as.numeric(v[present]))
      non_numeric <- mean(is.na(num))
      uniq <- length(unique(v[present])) / sum(present)
      if (non_numeric > 0.5 || uniq < 0.5) i <- i + 1 else break
    }
    anno_idx <- if (i > 2) 2:(i - 1) else integer(0)
    log_msg("Detected annotation rows: ",
            if (length(anno_idx)) paste(labels[anno_idx], collapse = ", ") else "none",
            ". Set input.annotation_rows to fix this if it is wrong.")
  } else {
    anno_idx <- match(annotation_rows, labels)
    if (anyNA(anno_idx))
      stop("Annotation row(s) not found: ",
           paste(annotation_rows[is.na(anno_idx)], collapse = ", "), call. = FALSE)
  }
  samples <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
  for (j in anno_idx) samples[[make.names(labels[j])]] <- trimws(as.character(unlist(raw[j, -1])))

  first_data <- max(c(1, anno_idx)) + 1
  block <- raw[first_data:nrow(raw), , drop = FALSE]
  m <- .to_numeric_matrix(block[-1])
  n_bad <- sum(is.na(m) & !is.na(as.matrix(block[-1])))
  if (n_bad) log_msg("NOTE: ", n_bad, " non-numeric data cells became NA.")
  colnames(m) <- ids
  out <- .finish(m, block[[1]], block[[1]])
  out$samples <- samples
  log_msg("Read ", nrow(m), " features x ", ncol(m), " columns (annotated_xlsx).")
  out
}
