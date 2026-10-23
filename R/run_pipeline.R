#' Run the full analysis
#'
#' Reads data, preprocesses, tests every configured comparison, runs Reactome
#' enrichment (if the species and packages allow) and AI interpretation (if
#' enabled and an API key is available), then writes results to
#' `output$dir`.
#'
#' Outputs: `DE_<comparison>.csv`, `DE_summary.csv`, `F_test.csv` (> 2 groups,
#' limma), `Reactome_ORA_*.csv` / `Reactome_GSEA_*.csv`, `PCA.png`,
#' `Volcano_*.png`, `AI_interpretation.md`, `config_used.yml`, `run_log.txt`,
#' `session_info.txt` and `results.rds`.
#'
#' @param config Path to a YAML file, or a configuration list.
#' @param verbose Print progress messages.
#' @return A list with `se`, `de`, `enrichment`, `ai` and `config`, invisibly.
#' @export
#' @examples
#' \donttest{
#' cfg_path <- simulate_dataset(tempfile("demo"), n_groups = 3)
#' res <- run_pipeline(cfg_path)
#' res$de$summary
#' }
run_pipeline <- function(config, verbose = TRUE) {
  .pkg_state$verbose <- verbose
  log_reset()
  cfg <- if (is.character(config)) read_config(config) else validate_config(config)
  out <- cfg$output$dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (file.access(out, 2) != 0) stop("Output folder is not writable: ", out, call. = FALSE)
  log_msg("sparcproteomics ", utils::packageVersion("sparcproteomics"), " | ",
          cfg$study$title, " | species: ", cfg$study$species)

  se <- build_dataset(cfg)
  se_unfiltered <- se            # detection testing needs proteins the filter removes
  pp <- cfg$preprocessing
  se <- pp_filter(se, pp$filter$method, pp$filter$min_valid_fraction)
  ref <- pp$normalization_reference
  if (!is.null(ref)) {
    if (!file.exists(ref)) stop("Normalization reference file not found: ", ref, call. = FALSE)
    ref <- trimws(readLines(ref, warn = FALSE))
    ref <- ref[nzchar(ref) & !startsWith(ref, "#")]
  }
  se <- pp_normalize(se, pp$normalization, reference = ref)
  im <- pp$imputation
  se <- pp_impute(se, im$method, shift = im$downshift_shift, width = im$downshift_width,
                  mnar_group_fraction = im$mnar_group_fraction, k = im$knn_k, seed = im$seed)

  comps <- build_comparisons(levels(se$group), cfg$design$comparisons, levels(se$group)[1])
  log_msg("Comparisons: ", paste(comps$name, collapse = ", "))
  st <- cfg$statistics
  de <- run_de(se, comps, method = st$method, covariates = cfg$design$covariates,
               block = cfg$design$block, min_valid_per_group = st$min_valid_per_group,
               p_adjust = st$p_adjust, adj_p_cutoff = st$adj_p_cutoff,
               log2fc_cutoff = st$log2fc_cutoff, trend = st$limma_trend,
               robust = st$limma_robust, global_f_test = st$global_f_test,
               significance = st$significance, treat_fc = st$treat_fc,
               use_imputed = im$use_for_de)

  .check_centering(de, pp$normalization)
  for (nm in names(de$tables))
    utils::write.csv(de$tables[[nm]], file.path(out, paste0("DE_", nm, ".csv")), row.names = FALSE)
  utils::write.csv(de$summary, file.path(out, "DE_summary.csv"), row.names = FALSE)
  if (!is.null(de$f_test))
    utils::write.csv(de$f_test, file.path(out, "F_test.csv"), row.names = FALSE)

  detection <- NULL
  if (isTRUE(cfg$detection$enabled)) {
    detection <- run_detection(se_unfiltered, comps, adj_p_cutoff = st$adj_p_cutoff,
                               min_delta = cfg$detection$min_delta,
                               on_off_min_fraction = cfg$detection$on_off_min_fraction)
    for (nm in names(detection$tables))
      utils::write.csv(detection$tables[[nm]], file.path(out, paste0("Detection_", nm, ".csv")),
                       row.names = FALSE)
    utils::write.csv(detection$summary, file.path(out, "Detection_summary.csv"), row.names = FALSE)
  }

  enrichment <- NULL
  if (isTRUE(cfg$enrichment$enabled)) {
    enrichment <- lapply(de$tables, run_enrichment, species = cfg$study$species,
                         min_gs_size = cfg$enrichment$min_gs_size,
                         max_gs_size = cfg$enrichment$max_gs_size,
                         p_cutoff = cfg$enrichment$p_cutoff)
    for (nm in names(enrichment)) {
      log_msg("Reactome ", nm, ": ", enrichment[[nm]]$status)
      if (nrow(enrichment[[nm]]$ora))
        utils::write.csv(enrichment[[nm]]$ora, file.path(out, paste0("Reactome_ORA_", nm, ".csv")),
                         row.names = FALSE)
      if (nrow(enrichment[[nm]]$gsea))
        utils::write.csv(enrichment[[nm]]$gsea, file.path(out, paste0("Reactome_GSEA_", nm, ".csv")),
                         row.names = FALSE)
    }
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    pca_assay <- if (im$method == "none") "norm" else "imputed"
    tryCatch(ggplot2::ggsave(file.path(out, "PCA.png"), plot_pca(se, pca_assay),
                             width = 6, height = 4.5, dpi = 300),
             error = function(e) log_msg("PCA plot skipped: ", conditionMessage(e)))
    for (nm in names(de$tables))
      ggplot2::ggsave(file.path(out, paste0("Volcano_", nm, ".png")),
                      plot_volcano(de$tables[[nm]], gsub("_", " ", nm), st$adj_p_cutoff,
                                   if (st$significance == "treat") log2(st$treat_fc) else st$log2fc_cutoff,
                                   subtitle = sprintf("%d significant (%s)", sum(de$tables[[nm]]$significant),
                                                      significance_label(de$settings))),
                      width = 6, height = 5, dpi = 300)
  }

  ai <- ai_interpret(de, se, cfg, enrichment, detection = detection)
  if (!is.null(ai)) writeLines(ai, file.path(out, "AI_interpretation.md"))

  write_config(cfg, file.path(out, "config_used.yml"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "session_info.txt"))
  res <- list(se = se, de = de, detection = detection, enrichment = enrichment, ai = ai,
              config = cfg)
  saveRDS(res, file.path(out, "results.rds"))
  log_msg("Done. Results in ", normalizePath(out))
  writeLines(get_run_log(), file.path(out, "run_log.txt"))
  invisible(res)
}

utils::globalVariables(".data")
