#' Build the analysis dataset from a configuration
#'
#' Reads the protein table and sample sheet, matches samples, resolves repeated
#' column names, checks the annotations the design needs, and returns a
#' `SummarizedExperiment`. Matrix and sample annotations are subset together so
#' they cannot drift out of alignment.
#'
#' Checks that stop the run (rather than silently producing a wrong design):
#' the group column or control level is missing; a covariate or block column
#' is missing; repeated sample names with `duplicate_columns = "error"`;
#' fewer than two groups remain.
#'
#' @param cfg Configuration from [read_config()] or [default_config()].
#' @return A `SummarizedExperiment` with assays `intensity` (as read) and
#'   `log2`, `colData` holding the sample sheet plus a `group` factor whose
#'   first level is the control, and `rowData` with `protein_id` and `gene`.
#' @export
build_dataset <- function(cfg) {
  cfg <- validate_config(cfg)
  inp <- cfg$input
  tab <- if (inp$format %in% .precursor_formats) .read_and_summarize(inp) else
    read_protein_table(inp$matrix, inp$format, intensity = inp$intensity,
                       id_column = inp$id_column, gene_column = inp$gene_column,
                       annotation_rows = inp$annotation_rows)

  if (inp$format == "annotated_xlsx") {
    samples <- tab$samples
    sample_col <- "sample_id"
  } else {
    samples <- read_sample_sheet(inp$sample_sheet)
    sample_col <- inp$sample_column
    if (!sample_col %in% names(samples))
      stop("Sample-ID column '", sample_col, "' not in the sample sheet. Columns: ",
           paste(names(samples), collapse = ", "), call. = FALSE)
    if (anyDuplicated(samples[[sample_col]]))
      stop("Sample sheet has repeated sample IDs: ",
           paste(unique(samples[[sample_col]][duplicated(samples[[sample_col]])]),
                 collapse = ", "), call. = FALSE)
  }

  m <- tab$intensities
  design_cols <- unique(c(cfg$design$group_column, unlist(cfg$design$covariates),
                          cfg$design$block))

  # repeated sample names (mostly annotated workbooks with several metrics per run)
  if (anyDuplicated(colnames(m))) {
    if (inp$format == "annotated_xlsx") {
      res <- resolve_duplicate_columns(m, samples, inp$duplicate_columns,
                                       intersect(design_cols, names(samples)))
      m <- res$matrix; samples <- res$samples
      if (identical(inp$duplicate_columns, "log2") && !isTRUE(inp$is_log2)) {
        log_msg("input.is_log2 set to TRUE: the retained columns are already log2.")
        inp$is_log2 <- TRUE
      }
    } else {
      stop("The protein table has repeated sample column names: ",
           paste(unique(colnames(m)[duplicated(colnames(m))]), collapse = ", "),
           call. = FALSE)
    }
  }

  # match sample sheet rows to matrix columns
  idx <- match_samples(colnames(m), samples[[sample_col]])
  unmatched_cols <- colnames(m)[is.na(idx)]
  if (length(unmatched_cols))
    log_msg("WARNING: ", length(unmatched_cols), " matrix column(s) not in the sample ",
            "sheet and dropped: ", paste(utils::head(unmatched_cols, 5), collapse = ", "),
            if (length(unmatched_cols) > 5) ", ..." else "")
  missing_rows <- setdiff(seq_len(nrow(samples)), idx)
  if (length(missing_rows))
    log_msg("WARNING: ", length(missing_rows), " sample-sheet row(s) have no matrix ",
            "column: ", paste(utils::head(samples[[sample_col]][missing_rows], 5),
                              collapse = ", "))
  keep <- !is.na(idx)
  if (sum(keep) < 2) stop("Fewer than two samples matched between the matrix and ",
                          "the sample sheet. Check input.sample_column.", call. = FALSE)
  m <- m[, keep, drop = FALSE]
  samples <- samples[idx[keep], , drop = FALSE]
  colnames(m) <- samples[[sample_col]]
  rownames(samples) <- samples[[sample_col]]

  samples <- .check_annotations(samples, cfg, design_cols)
  m <- m[, rownames(samples), drop = FALSE]

  # values
  if (isTRUE(inp$zero_is_missing)) {
    nz <- sum(m == 0, na.rm = TRUE)
    if (nz) log_msg("Converted ", nz, " zero values to NA (not detected).")
    m[!is.na(m) & m == 0] <- NA
  }
  if (!isTRUE(inp$is_log2) && any(m < 0, na.rm = TRUE))
    stop("Negative intensities found but input.is_log2 is FALSE. If the data are ",
         "already log-transformed, set input.is_log2: true.", call. = FALSE)
  if (!isTRUE(inp$is_log2) && stats::median(m, na.rm = TRUE) < 40)
    log_msg("WARNING: median intensity is ", signif(stats::median(m, na.rm = TRUE), 3),
            ", which looks log-scaled, but input.is_log2 is FALSE. Logging twice ",
            "compresses every fold change.")
  lg <- if (isTRUE(inp$is_log2)) m else log2(m)
  lg[!is.finite(lg)] <- NA

  feats <- tab$features[rownames(m), , drop = FALSE]
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = m, log2 = lg),
    rowData = S4Vectors::DataFrame(feats),
    colData = S4Vectors::DataFrame(samples, check.names = FALSE),
    metadata = list(config = cfg))
  log_msg("Dataset: ", nrow(se), " proteins x ", ncol(se), " samples; groups: ",
          paste(sprintf("%s (n=%d)", levels(se$group), as.integer(table(se$group))),
                collapse = ", "), "; control = '", levels(se$group)[1], "'.")
  se
}

.check_annotations <- function(samples, cfg, design_cols) {
  gcol <- cfg$design$group_column
  miss <- setdiff(design_cols, names(samples))
  if (length(miss))
    stop("Design column(s) not in the sample annotations: ",
         paste(miss, collapse = ", "), ". Available: ",
         paste(names(samples), collapse = ", "), call. = FALSE)

  blank <- function(x) is.na(x) | !nzchar(trimws(x)) |
    tolower(trimws(x)) %in% c("na", "n/a", "nan", "null", "none", "-", "?")
  bad <- Reduce(`|`, lapply(design_cols, function(k) blank(samples[[k]])))
  if (any(bad)) {
    msg <- paste0(sum(bad), " sample(s) have blank design annotations: ",
                  paste(utils::head(rownames(samples)[bad], 5), collapse = ", "))
    if (identical(cfg$design$incomplete_annotation, "error")) stop(msg, call. = FALSE)
    log_msg("WARNING: ", msg, ". Dropped (design.incomplete_annotation = drop).")
    samples <- samples[!bad, , drop = FALSE]
  }

  g <- samples[[gcol]]
  ctrl <- cfg$design$control
  if (is.null(ctrl)) {
    ctrl <- sort(unique(g))[1]
    log_msg("WARNING: design.control not set; using '", ctrl, "' as the control.")
  }
  if (!ctrl %in% g)
    stop("Control group '", ctrl, "' not found in column '", gcol, "'. Levels: ",
         paste(unique(g), collapse = ", "), call. = FALSE)

  n <- table(g)
  small <- names(n)[n < cfg$design$min_group_size]
  if (ctrl %in% small)
    stop("Control group '", ctrl, "' has only ", n[[ctrl]], " sample(s).", call. = FALSE)
  if (length(small)) {
    log_msg("WARNING: group(s) below design.min_group_size (",
            cfg$design$min_group_size, ") dropped: ", paste(small, collapse = ", "))
    samples <- samples[!g %in% small, , drop = FALSE]
    g <- samples[[gcol]]
  }
  lv <- c(ctrl, setdiff(unique(g), ctrl))
  if (length(lv) < 2) stop("Fewer than two groups remain.", call. = FALSE)
  samples$group <- factor(g, levels = lv)

  for (k in unlist(cfg$design$covariates)) {
    v <- suppressWarnings(as.numeric(samples[[k]]))
    samples[[k]] <- if (!anyNA(v) && length(unique(v)) > 5) v else factor(samples[[k]])
  }
  if (!is.null(cfg$design$block)) samples[[cfg$design$block]] <- factor(samples[[cfg$design$block]])
  samples
}

#' Match matrix column names to sample-sheet IDs
#'
#' Tries an exact match first, then a relaxed match that ignores folders, file
#' extensions (`.raw`, `.d`, `.wiff`, `.mzML` ...), Spectronaut `[n]` prefixes
#' and case. A relaxed key that matches more than one sample is not used.
#'
#' @param cols Matrix column names.
#' @param ids Sample-sheet IDs.
#' @return Integer vector: for each column, the matching row of the sample
#'   sheet, or `NA`.
#' @export
#' @examples
#' match_samples(c("D:/runs/S01.raw", "[2] S02.d"), c("S02", "S01"))
match_samples <- function(cols, ids) {
  idx <- match(cols, ids)
  key <- function(x) {
    x <- basename(gsub("\\\\", "/", x))
    x <- sub("^\\[[0-9]+\\]\\s*", "", x)
    x <- sub("\\.PG\\.Quantity$", "", x)
    x <- sub("\\.(raw|d|wiff2?|mzml|mzxml|dia|htrms)$", "", x, ignore.case = TRUE)
    tolower(trimws(x))
  }
  kc <- key(cols); ki <- key(ids)
  ki[duplicated(ki) | duplicated(ki, fromLast = TRUE)] <- NA
  relaxed <- match(kc, ki)
  idx[is.na(idx)] <- relaxed[is.na(idx)]
  idx[duplicated(idx) & !is.na(idx)] <- NA
  idx
}

#' Resolve repeated sample columns
#'
#' Repeated sample names must be resolved before analysis: left alone, each
#' copy is counted as an independent replicate and every p-value is too small.
#'
#' @param m Numeric matrix with repeated column names.
#' @param samples Sample annotations, one row per column of `m`, with
#'   `sample_id`.
#' @param mode `"error"` stops; `"first"` keeps the first column per name;
#'   `"average"` averages columns sharing a name; `"log2"` keeps, per name, the
#'   column with the smallest median (the log2 column when an export holds
#'   linear quantities and their log2 side by side); `"suffix"` makes names
#'   unique and treats copies as independent (diagnostic use only).
#' @param check_cols Annotation columns that must agree within a name before
#'   columns are collapsed.
#' @return A list with `matrix` and `samples`.
#' @export
resolve_duplicate_columns <- function(m, samples, mode = "error", check_cols = character(0)) {
  ids <- colnames(m)
  grp <- factor(ids, levels = unique(ids))
  k <- max(table(grp))
  log_msg("Repeated sample names: ", sum(table(grp) > 1), " names appear up to ",
          k, " times.")
  if (mode == "error")
    stop("Repeated sample names in the header. Set input.duplicate_columns to ",
         "first, average, log2 or suffix once you know why they repeat.", call. = FALSE)
  if (mode == "suffix") {
    log_msg("WARNING: duplicate_columns = suffix treats repeated columns as ",
            "independent samples; p-values will be anticonservative.")
    new <- make.unique(ids)
    colnames(m) <- new; samples$sample_id <- new
    return(list(matrix = m, samples = samples))
  }
  for (cc in check_cols) {
    clash <- tapply(samples[[cc]], grp, function(v) length(unique(v)) > 1)
    if (any(clash))
      stop("Columns sharing a sample name disagree on '", cc, "' (e.g. ",
           names(clash)[clash][1], "); they are not copies of one sample.", call. = FALSE)
  }
  if (mode == "average") {
    out <- vapply(levels(grp), function(g)
      rowMeans(m[, grp == g, drop = FALSE], na.rm = TRUE), numeric(nrow(m)))
    out[is.nan(out)] <- NA
    rownames(out) <- rownames(m)
    return(list(matrix = out, samples = samples[!duplicated(ids), , drop = FALSE]))
  }
  sel <- if (mode == "first") which(!duplicated(ids)) else {
    med <- apply(m, 2, function(x) stats::median(x[is.finite(x)], na.rm = TRUE))
    vapply(split(seq_along(ids), grp), function(ix)
      if (all(is.na(med[ix]))) ix[1] else ix[which.min(med[ix])], integer(1))
  }
  sel <- sort(sel)
  log_msg("Kept ", length(sel), " of ", length(ids), " columns (mode = ", mode, ").")
  list(matrix = m[, sel, drop = FALSE], samples = samples[sel, , drop = FALSE])
}
