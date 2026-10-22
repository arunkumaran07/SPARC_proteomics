#' Simulate a demo dataset
#'
#' Writes a simulated protein table in DIA-NN, Spectronaut or MaxQuant format,
#' a matching sample sheet, the true differential proteins, and a ready-to-run
#' YAML configuration. Useful for trying options and for tests; contains no
#' real data.
#'
#' Log2 abundances are drawn per protein, each non-control group shifts a set
#' of proteins by `effect` log2 units (half up, half down), and values go
#' missing preferentially at low abundance (detection limit) plus a small
#' random fraction.
#'
#' @param dir Output folder (created).
#' @param n_groups Number of groups including the control.
#' @param n_per_group Samples per group.
#' @param n_proteins Number of proteins.
#' @param frac_de Fraction of proteins changed in each non-control group.
#' @param effect Size of the change in log2 units.
#' @param format `"diann"`, `"spectronaut"` or `"maxquant"`.
#' @param species Written into the configuration.
#' @param seed Random seed.
#' @return Path to the configuration file.
#' @export
#' @examples
#' cfg <- simulate_dataset(tempfile("demo"), n_groups = 3)
#' readLines(cfg)[1:10]
simulate_dataset <- function(dir, n_groups = 2, n_per_group = 4, n_proteins = 1500,
                             frac_de = 0.05, effect = 1.5,
                             format = c("diann", "spectronaut", "maxquant"),
                             species = "human", seed = 1) {
  format <- match.arg(format)
  stopifnot(n_groups >= 2, n_per_group >= 2)
  set.seed(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  groups <- c("Control", paste0("Treatment", LETTERS[seq_len(n_groups - 1)]))
  ss <- data.frame(sample_id = sprintf("S%02d", seq_len(n_groups * n_per_group)),
                   group = rep(groups, each = n_per_group),
                   sex = rep(c("F", "M"), length.out = n_groups * n_per_group),
                   stringsAsFactors = FALSE)
  base <- stats::rnorm(n_proteins, 22, 2)
  load <- stats::rnorm(nrow(ss), 0, 0.2)
  x <- outer(base, rep(1, nrow(ss))) + matrix(stats::rnorm(n_proteins * nrow(ss), 0, 0.35),
                                              n_proteins) + rep(load, each = n_proteins)
  ids <- sprintf("P%05d", seq_len(n_proteins))
  genes <- sprintf("GENE%04d", seq_len(n_proteins))
  truth <- list()
  n_de <- round(frac_de * n_proteins)
  for (g in groups[-1]) {
    if (n_de == 0) next
    hit <- sample(n_proteins, n_de)
    eff <- rep(c(effect, -effect), length.out = n_de)
    x[hit, ss$group == g] <- x[hit, ss$group == g] + eff
    truth[[g]] <- data.frame(group = g, protein_id = ids[hit], true_log2fc = eff)
  }
  p_miss <- stats::plogis(-(x - 18.5) * 1.5) + 0.02
  x[matrix(stats::runif(length(x)), nrow(x)) < p_miss] <- NA
  lin <- 2^x

  path <- switch(format,
    diann = {
      df <- data.frame(Protein.Group = ids, Protein.Ids = ids,
                       Protein.Names = paste0(genes, "_HUMAN"), Genes = genes,
                       First.Protein.Description = paste("Simulated protein", seq_len(n_proteins)),
                       check.names = FALSE)
      m <- lin; colnames(m) <- paste0("D:/raw/", ss$sample_id, ".raw")
      file.path(dir, "report.pg_matrix.tsv") -> f
      utils::write.table(cbind(df, m), f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      f
    },
    spectronaut = {
      df <- data.frame(PG.ProteinGroups = ids, PG.Genes = genes,
                       PG.ProteinDescriptions = paste("Simulated protein", seq_len(n_proteins)),
                       check.names = FALSE)
      m <- lin; m[is.na(m)] <- NaN
      colnames(m) <- sprintf("[%d] %s.raw.PG.Quantity", seq_len(nrow(ss)), ss$sample_id)
      f <- file.path(dir, "spectronaut_protein_report.tsv")
      utils::write.table(cbind(df, m), f, sep = "\t", quote = FALSE, row.names = FALSE)
      f
    },
    maxquant = {
      df <- data.frame(`Protein IDs` = ids, `Gene names` = genes,
                       `Protein names` = paste("Simulated protein", seq_len(n_proteins)),
                       Reverse = "", `Potential contaminant` = "", check.names = FALSE)
      df$Reverse[1:5] <- "+"
      m <- lin; m[is.na(m)] <- 0
      colnames(m) <- paste("LFQ intensity", ss$sample_id)
      f <- file.path(dir, "proteinGroups.txt")
      utils::write.table(cbind(df, m), f, sep = "\t", quote = FALSE, row.names = FALSE)
      f
    })
  utils::write.csv(ss, file.path(dir, "samples.csv"), row.names = FALSE)
  truth <- if (length(truth)) do.call(rbind, truth) else
    data.frame(group = character(0), protein_id = character(0), true_log2fc = numeric(0))
  utils::write.csv(truth, file.path(dir, "truth.csv"), row.names = FALSE)

  cfg <- default_config()
  cfg$study$title <- sprintf("Simulated %d-group study", n_groups)
  cfg$study$species <- species
  cfg$input$matrix <- basename(path)
  cfg$input$format <- format
  cfg$input$sample_sheet <- "samples.csv"
  cfg$design$group_column <- "group"
  cfg$design$control <- "Control"
  cfg$output$dir <- "results"
  cfg$ai$enabled <- FALSE
  cp <- file.path(dir, "config.yml")
  write_config(cfg, cp)
  cp
}
