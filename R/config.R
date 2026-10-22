#' Default analysis configuration
#'
#' Every option the pipeline understands, with its default. A YAML file read by
#' [read_config()] only needs to contain the values that differ from these.
#'
#' The main choices:
#' \describe{
#'   \item{`study$species`}{See [list_species()]. Controls Reactome enrichment.}
#'   \item{`design$group_column`, `design$control`}{Which sample-sheet column
#'     defines the groups, and which level is the control (reference).}
#'   \item{`design$comparisons`}{`"all_vs_control"`, `"pairwise"`, or a list of
#'     explicit comparisons written `"A - B"` or `"A vs B"` (A is the numerator,
#'     so positive log2FC means higher in A).}
#'   \item{`preprocessing$imputation$method`}{`"none"`, `"mixed"`, `"knn"`,
#'     `"downshift"` or `"mindet"`; see [pp_impute()].}
#'   \item{`statistics$method`}{`"limma"`, `"welch"` or `"student"`; see
#'     [run_de()].}
#'   \item{`ai$enabled`}{`"auto"` (use AI only if an API key is set), `TRUE`
#'     or `FALSE`; see [ai_interpret()].}
#' }
#'
#' @return A nested list.
#' @export
#' @examples
#' str(default_config(), max.level = 1)
default_config <- function() {
  list(
    study = list(
      title       = "Proteomics analysis",
      description = "",      # free text; used in the AI prompt and the report
      species     = "human", # see list_species()
      tissue      = ""       # free text, e.g. "plasma", "skeletal muscle"
    ),
    input = list(
      matrix            = NULL,      # path to the protein-group table
      format            = "diann",   # diann | spectronaut | maxquant | generic | annotated_xlsx
      sample_sheet      = NULL,      # CSV/TSV/XLSX; not needed for annotated_xlsx
      sample_column     = "sample_id",
      intensity         = "auto",    # maxquant: "lfq" | "intensity"; auto prefers LFQ
      id_column         = NULL,      # generic format only
      gene_column       = NULL,      # generic format only
      is_log2           = FALSE,
      zero_is_missing   = TRUE,
      annotation_rows   = NULL,      # annotated_xlsx: row labels holding sample annotations
      duplicate_columns = "error"    # error | first | average | suffix | log2
    ),
    output = list(dir = "results"),
    design = list(
      group_column          = "group",
      control               = NULL,       # reference level of group_column
      comparisons           = "all_vs_control",
      covariates            = list(),     # sample-sheet columns to adjust for
      block                 = NULL,       # e.g. animal/patient ID for repeated measures
      incomplete_annotation = "drop",     # drop | error
      min_group_size        = 2
    ),
    preprocessing = list(
      filter = list(
        method             = "any_group", # any_group | all_groups | overall | none
        min_valid_fraction = 0.7
      ),
      normalization = "median",           # none | median | quantile | cyclicloess
      imputation = list(
        method              = "mixed",    # none | mixed | knn | downshift | mindet
        downshift_shift     = 1.8,        # SDs below the sample mean
        downshift_width     = 0.3,        # width as a fraction of the sample SD
        mnar_group_fraction = 0.6,        # mixed: >= this missing in a group -> MNAR
        knn_k               = 10,
        seed                = 42,
        use_for_de          = FALSE       # TRUE tests the imputed matrix (not recommended)
      )
    ),
    statistics = list(
      method              = "limma",      # limma | welch | student
      limma_trend         = TRUE,
      limma_robust        = TRUE,
      p_adjust            = "BH",         # any method accepted by stats::p.adjust
      adj_p_cutoff        = 0.05,
      log2fc_cutoff       = 0.585,        # 1.5-fold
      min_valid_per_group = 2,            # per comparison, both groups need this many values
      global_f_test       = TRUE          # limma only, when there are > 2 groups
    ),
    enrichment = list(
      enabled      = TRUE,
      min_gs_size  = 10,
      max_gs_size  = 500,
      p_cutoff     = 0.05
    ),
    ai = list(
      enabled     = "auto",               # auto | TRUE | FALSE
      key_env     = "ANTHROPIC_API_KEY",  # environment variable holding the key
      model       = "claude-opus-5-5",
      max_tokens  = 16000,
      web_search  = FALSE,
      search_tool = "web_search_20260318",
      max_uses    = 20,
      max_turns   = 8,
      api_url     = "https://api.anthropic.com/v1/messages",
      api_version = "2023-06-01"
    )
  )
}

.choices <- list(
  format        = c("diann", "spectronaut", "maxquant", "generic", "annotated_xlsx"),
  duplicate     = c("error", "first", "average", "suffix", "log2"),
  incomplete    = c("drop", "error"),
  filter        = c("any_group", "all_groups", "overall", "none"),
  normalization = c("none", "median", "quantile", "cyclicloess"),
  imputation    = c("none", "mixed", "knn", "downshift", "mindet"),
  stat_method   = c("limma", "welch", "student")
)

#' Read a YAML configuration
#'
#' Values in the file override [default_config()]; anything not given keeps its
#' default. Relative file paths are resolved against the YAML file's folder.
#'
#' @param path Path to a YAML file.
#' @return A validated configuration list.
#' @export
read_config <- function(path) {
  if (!file.exists(path)) stop("Config file not found: ", path, call. = FALSE)
  user <- yaml::read_yaml(path)
  cfg <- utils::modifyList(default_config(), user %||% list())
  base <- dirname(normalizePath(path))
  for (k in c("matrix", "sample_sheet"))
    if (!is.null(cfg$input[[k]])) cfg$input[[k]] <- .resolve_path(cfg$input[[k]], base)
  cfg$output$dir <- .resolve_path(cfg$output$dir, base)
  validate_config(cfg)
}

#' Write a configuration to YAML
#'
#' @param cfg A configuration list.
#' @param path Output path.
#' @return `path`, invisibly.
#' @export
write_config <- function(cfg, path) {
  yaml::write_yaml(cfg, path)
  invisible(path)
}

#' Check a configuration
#'
#' Checks option values against the allowed choices. Checks that depend on the
#' data (does the control level exist?) happen in [build_dataset()].
#'
#' @param cfg A configuration list.
#' @return `cfg`, with species normalised, invisibly; stops on an invalid value.
#' @export
validate_config <- function(cfg) {
  chk <- function(value, allowed, name) {
    if (length(value) != 1 || !as.character(value) %in% allowed)
      stop("Config '", name, "' is '", paste(value, collapse = ","),
           "'; must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  cfg$study$species <- species_info(cfg$study$species)$species
  chk(cfg$input$format, .choices$format, "input.format")
  chk(cfg$input$duplicate_columns, .choices$duplicate, "input.duplicate_columns")
  chk(cfg$design$incomplete_annotation, .choices$incomplete, "design.incomplete_annotation")
  chk(cfg$preprocessing$filter$method, .choices$filter, "preprocessing.filter.method")
  chk(cfg$preprocessing$normalization, .choices$normalization, "preprocessing.normalization")
  chk(cfg$preprocessing$imputation$method, .choices$imputation,
      "preprocessing.imputation.method")
  chk(cfg$statistics$method, .choices$stat_method, "statistics.method")
  chk(cfg$statistics$p_adjust, stats::p.adjust.methods, "statistics.p_adjust")
  chk(tolower(as.character(cfg$ai$enabled)), c("auto", "true", "false"), "ai.enabled")
  f <- cfg$preprocessing$filter$min_valid_fraction
  if (!is.numeric(f) || f < 0 || f > 1)
    stop("preprocessing.filter.min_valid_fraction must be between 0 and 1.", call. = FALSE)
  if (cfg$statistics$min_valid_per_group < 2)
    stop("statistics.min_valid_per_group must be at least 2.", call. = FALSE)
  if (is.null(cfg$input$matrix))
    stop("input.matrix (path to the protein table) is required.", call. = FALSE)
  if (cfg$input$format != "annotated_xlsx" && is.null(cfg$input$sample_sheet))
    stop("input.sample_sheet is required for format '", cfg$input$format, "'.",
         call. = FALSE)
  invisible(cfg)
}

#' Create a configuration interactively
#'
#' Asks, in the R console, for the species, input files and format, the group
#' column, the control group, the comparisons, the imputation and statistical
#' methods, and whether to use AI interpretation. It reads the sample sheet so
#' groups are chosen from the levels actually present. The result is written
#' to YAML so the same analysis can be re-run non-interactively with
#' [run_pipeline()].
#'
#' @param path Where to write the YAML file.
#' @param matrix,sample_sheet Optional file paths; asked for if `NULL`.
#' @return The configuration list, invisibly.
#' @export
create_config <- function(path = "analysis_config.yml", matrix = NULL,
                          sample_sheet = NULL) {
  if (!interactive())
    stop("create_config() needs an interactive R session. Edit a YAML file ",
         "and use read_config() instead.", call. = FALSE)
  cfg <- default_config()
  pick <- function(choices, title, multiple = FALSE, labels = choices) {
    sel <- utils::select.list(labels, title = title, multiple = multiple,
                              graphics = FALSE)
    if (!length(sel) || identical(sel, "")) stop("Cancelled.", call. = FALSE)
    choices[match(sel, labels)]
  }
  ask <- function(prompt, default = "") {
    ans <- readline(sprintf("%s%s: ", prompt,
                            if (nzchar(default)) paste0(" [", default, "]") else ""))
    if (nzchar(ans)) ans else default
  }

  cfg$study$title <- ask("Study title", cfg$study$title)
  cfg$study$description <- ask("One-paragraph study description (optional)")
  sp <- list_species()
  cfg$study$species <- pick(sp$species, "Species",
                            labels = paste0(sp$species, " (", sp$label, ")"))
  cfg$study$tissue <- ask("Tissue / sample type (optional)")

  cfg$input$format <- pick(.choices$format, "Input format")
  cfg$input$matrix <- matrix %||% ask("Path to the protein table")
  if (cfg$input$format != "annotated_xlsx") {
    cfg$input$sample_sheet <- sample_sheet %||% ask("Path to the sample sheet")
    ss <- read_sample_sheet(cfg$input$sample_sheet)
    cfg$input$sample_column <- pick(names(ss), "Column holding sample IDs")
  } else {
    ss <- read_annotated_xlsx(cfg$input$matrix)$samples
    cfg$input$sample_column <- "sample_id"
  }
  cfg$input$is_log2 <- pick(c(FALSE, TRUE), "Are the intensities already log2?",
                            labels = c("No (linear)", "Yes (log2)")) == "TRUE"

  design_cols <- setdiff(names(ss), cfg$input$sample_column)
  cfg$design$group_column <- pick(design_cols, "Column defining the sample groups")
  lv <- unique(as.character(ss[[cfg$design$group_column]]))
  lv <- lv[!is.na(lv) & nzchar(lv)]
  cat("Groups found:", paste(lv, collapse = ", "), "\n")
  cfg$design$control <- pick(lv, "Control (reference) group")

  mode <- pick(c("all_vs_control", "pairwise", "custom"), "Comparisons",
               labels = c("Every group vs the control", "All pairwise comparisons",
                          "Choose specific comparisons"))
  if (mode == "custom") {
    pairs <- utils::combn(lv, 2, simplify = FALSE)
    lab <- unlist(lapply(pairs, function(p) c(paste(p[2], "vs", p[1]),
                                              paste(p[1], "vs", p[2]))))
    cfg$design$comparisons <- as.list(pick(lab, "Select comparisons (A vs B = A relative to B)",
                                           multiple = TRUE))
  } else cfg$design$comparisons <- mode

  others <- setdiff(design_cols, cfg$design$group_column)
  if (length(others)) {
    cv <- utils::select.list(others, multiple = TRUE, graphics = FALSE,
                             title = "Covariates to adjust for (none = press Enter)")
    cfg$design$covariates <- as.list(cv)
    bl <- utils::select.list(c("(none)", others), graphics = FALSE,
                             title = "Blocking variable for repeated measures (e.g. animal ID)")
    cfg$design$block <- if (bl %in% c("", "(none)")) NULL else bl
  }

  cfg$preprocessing$normalization <- pick(.choices$normalization, "Normalization",
    labels = c("none (already normalized)", "median centering",
               "quantile", "cyclic loess"))
  cfg$preprocessing$imputation$method <- pick(.choices$imputation,
    "Imputation (used for PCA / clustering; tests use un-imputed data by default)",
    labels = c("none", "mixed: kNN for sporadic + down-shift for group-wise missing",
               "kNN only", "down-shift (Perseus-style, per sample)",
               "minimum per sample"))
  cfg$statistics$method <- pick(.choices$stat_method, "Statistical test",
    labels = c("limma moderated t-test (recommended)", "Welch t-test",
               "Student t-test"))

  has_key <- nzchar(Sys.getenv(cfg$ai$key_env))
  cat(if (has_key) "An Anthropic API key was found in " else "No API key found in ",
      cfg$ai$key_env, ".\n", sep = "")
  cfg$ai$enabled <- pick(c("auto", "TRUE", "FALSE"), "AI interpretation",
    labels = c("Auto: use it only if an API key is set",
               "Yes (requires the API key)", "No"))
  if (cfg$ai$enabled != "auto") cfg$ai$enabled <- as.logical(cfg$ai$enabled)

  cfg$output$dir <- ask("Output folder", cfg$output$dir)
  validate_config(cfg)
  write_config(cfg, path)
  message("Configuration written to ", normalizePath(path, mustWork = FALSE),
          "\nRun it with: run_pipeline(\"", path, "\")")
  invisible(cfg)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.resolve_path <- function(p, base) {
  if (is.null(p) || !nzchar(p)) return(p)
  if (grepl("^([A-Za-z]:|/|~)", p)) return(path.expand(p))
  file.path(base, p)
}
