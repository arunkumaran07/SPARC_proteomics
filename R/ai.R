#' Is AI interpretation available?
#'
#' AI interpretation is optional. With `ai$enabled = "auto"` it runs only when
#' the API key environment variable (default `ANTHROPIC_API_KEY`) is set and
#' \pkg{httr2} is installed. The key is read from the environment and never
#' written to any output.
#'
#' @param cfg Configuration list.
#' @return A list with `use` (logical) and `reason`.
#' @export
ai_status <- function(cfg) {
  en <- tolower(as.character(cfg$ai$enabled))
  if (en == "false") return(list(use = FALSE, reason = "disabled in config"))
  has_key <- nzchar(Sys.getenv(cfg$ai$key_env))
  has_pkg <- requireNamespace("httr2", quietly = TRUE)
  if (has_key && has_pkg) return(list(use = TRUE, reason = "API key found"))
  why <- c(if (!has_key) paste0("no API key in ", cfg$ai$key_env),
           if (!has_pkg) "package httr2 not installed")
  if (en == "true")
    log_msg("WARNING: ai.enabled is TRUE but ", paste(why, collapse = " and "),
            "; skipping AI interpretation.")
  list(use = FALSE, reason = paste(why, collapse = "; "))
}

#' Interpret results with Claude (optional)
#'
#' Sends only summary statistics (group sizes, comparison counts, top proteins
#' and pathways), never the raw data, and asks for a cautious interpretation.
#' The prompt is built from the configuration, so it fits any tissue, species
#' and design. Guard statements are added when results are null or
#' underpowered so the text cannot over-claim. The output should be reviewed
#' by a scientist before it goes into any report.
#'
#' @param de Result of [run_de()].
#' @param se The analysed `SummarizedExperiment`.
#' @param cfg Configuration list.
#' @param enrichment Optional named list of [run_enrichment()] results.
#' @param top_n Proteins listed per comparison.
#' @param detection Optional result of [run_detection()].
#' @return Markdown text, or `NULL` if AI is unavailable or the call failed.
#' @export
ai_interpret <- function(de, se, cfg, enrichment = NULL, top_n = 15, detection = NULL) {
  st <- ai_status(cfg)
  if (!st$use) { log_msg("AI interpretation skipped (", st$reason, ")."); return(NULL) }
  prompts <- .ai_prompts(de, se, cfg, enrichment, top_n, detection)
  txt <- .call_claude(prompts$system, prompts$user, cfg$ai)
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  paste0("> AI-generated interpretation (", cfg$ai$model, "). Review before use.\n\n", txt)
}

.ai_prompts <- function(de, se, cfg, enrichment, top_n, detection = NULL) {
  g <- se$group
  s <- de$summary
  tops <- vapply(names(de$tables), function(nm) {
    t <- de$tables[[nm]]
    t <- t[!is.na(t$adj.P.Val), , drop = FALSE][seq_len(min(top_n, sum(!is.na(t$adj.P.Val)))), ]
    paste0(nm, ":\n", paste(sprintf("  %s (%s): log2FC %.2f, adj.p %.2g%s",
      ifelse(is.na(t$gene), t$protein_id, t$gene), t$protein_id, t$logFC, t$adj.P.Val,
      ifelse(t$significant, " *", "")), collapse = "\n"))
  }, "")
  enr <- if (length(enrichment)) vapply(names(enrichment), function(nm) {
    e <- enrichment[[nm]]
    terms <- c(utils::head(e$ora$Description, 5), utils::head(e$gsea$Description, 5))
    paste0(nm, " (", e$status, "): ", if (length(terms)) paste(unique(terms), collapse = "; ") else "none")
  }, "") else "Not run."

  guards <- character(0)
  if (sum(s$significant) == 0) guards <- c(guards,
    "No protein passed the significance thresholds in any comparison. This is a null result: do not describe any protein as up- or down-regulated, do not build a mechanistic narrative, and discuss power and effect size instead.")
  if (min(table(g)) < 4) guards <- c(guards,
    "At least one group has fewer than 4 samples; emphasise limited power and that absence of significance is not evidence of no effect.")
  if (identical(de$settings$assay, "imputed")) guards <- c(guards,
    "Tests were run on imputed values; note that imputation can inflate significance.")

  system <- paste(
    "You are an expert proteomics and mass spectrometry scientist writing the",
    "interpretation section of an analysis report. Use only the numbers given.",
    "Never invent proteins, values, pathways or citations. Distinguish clearly",
    "between statistically significant results (marked *) and leading but",
    "non-significant ones. Write in plain scientific prose with three headed",
    "sections in Markdown: '## Summary', '## Interpretation', '## Limitations",
    "and next steps'.",
    if (isTRUE(cfg$ai$web_search)) "You may use web search to check protein functions; cite sources with DOI or PMID." else
      "Do not cite literature.")
  user <- paste0(
    "Study: ", cfg$study$title, "\n",
    if (nzchar(cfg$study$description)) paste0("Description: ", cfg$study$description, "\n"),
    "Species: ", species_info(cfg$study$species)$label, "\n",
    if (nzchar(cfg$study$tissue)) paste0("Sample type: ", cfg$study$tissue, "\n"),
    "Groups (n): ", paste(sprintf("%s (%d)", levels(g), as.integer(table(g))), collapse = ", "),
    "; control = ", levels(g)[1], "\n",
    "Proteins analysed: ", nrow(se), "\n",
    "Normalization: ", cfg$preprocessing$normalization,
    "; imputation (not used for testing unless stated): ", cfg$preprocessing$imputation$method, "\n",
    "Test: ", de$method, "; significant = ", significance_label(de$settings), "\n",
    if (length(de$settings$covariates)) paste0("Adjusted for: ", paste(de$settings$covariates, collapse = ", "), "\n"),
    "\nComparison summary:\n", paste(utils::capture.output(print(s, row.names = FALSE)), collapse = "\n"),
    if (!is.null(de$f_test)) paste0("\nGlobal F-test: ", sum(de$f_test$adj.P.Val < 0.05, na.rm = TRUE),
                                   " proteins differ across groups at BH < 0.05.\n") else "\n",
    "\nTop proteins per comparison (* = significant):\n", paste(tops, collapse = "\n\n"),
    "\n\nReactome enrichment:\n", paste(enr, collapse = "\n"),
    if (!is.null(detection)) paste0("\n\nPresence/absence (Fisher exact on detection; on/off = ",
      "descriptive, detected in one group only):\n", paste(utils::capture.output(
        print(detection$summary, row.names = FALSE)), collapse = "\n")) else "",
    if (length(guards)) paste0("\n\nCONSTRAINTS:\n- ", paste(guards, collapse = "\n- ")) else "")
  list(system = system, user = user)
}

.call_claude <- function(system, user, ai) {
  key <- Sys.getenv(ai$key_env)
  msgs <- list(list(role = "user", content = user))
  tools <- if (isTRUE(ai$web_search))
    list(list(type = ai$search_tool, name = "web_search", max_uses = ai$max_uses)) else NULL
  out <- character(0)
  for (turn in seq_len(ai$max_turns)) {
    body <- list(model = ai$model, max_tokens = ai$max_tokens, system = system,
                 messages = msgs)
    if (!is.null(tools)) body$tools <- tools
    resp <- tryCatch(
      httr2::request(ai$api_url) |>
        httr2::req_headers("x-api-key" = key, "anthropic-version" = ai$api_version,
                           "content-type" = "application/json") |>
        httr2::req_body_json(body, auto_unbox = TRUE) |>
        httr2::req_timeout(900) |>
        httr2::req_retry(max_tries = 4, backoff = function(i) min(60, 2^i)) |>
        httr2::req_perform() |>
        httr2::resp_body_json(),
      error = function(e) { log_msg("AI request failed: ", conditionMessage(e)); NULL })
    if (is.null(resp)) break
    txt <- paste(vapply(resp$content, function(b)
      if (identical(b$type, "text")) b$text else "", character(1)), collapse = "")
    out <- c(out, txt)
    log_msg("AI turn ", turn, ": stop_reason = ", resp$stop_reason)
    if (!identical(resp$stop_reason, "pause_turn")) break
    msgs <- c(msgs, list(list(role = "assistant", content = resp$content)))
  }
  trimws(paste(out, collapse = ""))
}
