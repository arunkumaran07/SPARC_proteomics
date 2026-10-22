# Regenerates the example datasets in inst/extdata/. Not part of the installed
# package (excluded in .Rbuildignore). Needs org.Hs.eg.db.
#
# The data are SIMULATED. Real human gene symbols and UniProt accessions are
# used so that annotation and Reactome enrichment behave as they would on real
# data, and two biologically coherent gene sets are changed on purpose:
#   TreatmentA: proteins annotated to electron transport chain (GO:0022900) up
#   TreatmentB: proteins annotated to cytoplasmic translation (GO:0002181) down
# plus a random set of proteins changed in each treatment. The planted changes
# are in truth.csv.

suppressPackageStartupMessages(library(org.Hs.eg.db))
set.seed(2026)
out <- file.path("inst", "extdata")

go_genes <- function(go) {
  g <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, go, "SYMBOL", "GOALL"))
  unique(g$SYMBOL)
}
pc <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, "protein-coding",
                                             c("SYMBOL", "UNIPROT"), "GENETYPE"))
pc <- pc[!is.na(pc$UNIPROT) & grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", pc$UNIPROT), ]
pc <- pc[!duplicated(pc$SYMBOL), ]

etc  <- intersect(go_genes("GO:0022900"), pc$SYMBOL)
tran <- intersect(go_genes("GO:0002181"), pc$SYMBOL)
etc  <- sample(etc, min(45, length(etc)))
tran <- sample(setdiff(tran, etc), min(45, length(tran)))
rest <- sample(setdiff(pc$SYMBOL, c(etc, tran)), 1800 - length(etc) - length(tran))
genes <- sample(c(etc, tran, rest))
acc <- pc$UNIPROT[match(genes, pc$SYMBOL)]
desc <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, genes, "GENENAME", "SYMBOL"))
n <- length(genes)

# design: 3 groups x 5, balanced sex
ss <- data.frame(sample_id = sprintf("HEK_%02d", 1:15),
                 group = rep(c("Control", "TreatmentA", "TreatmentB"), each = 5),
                 sex = rep(c("F", "M", "F", "M", "F"), 3),
                 batch = rep(c("B1", "B2"), length.out = 15),
                 stringsAsFactors = FALSE)
base <- rnorm(n, 21, 2.2)
x <- base + matrix(rnorm(n * 15, 0, 0.3), n) +
  rep(rnorm(15, 0, 0.25), each = n) +                 # loading differences
  rep(ifelse(ss$batch == "B2", 0.15, 0), each = n)    # mild batch effect

truth <- list()
spike <- function(idx, grp, eff) {
  x[idx, ss$group == grp] <<- x[idx, ss$group == grp] + eff
  truth[[length(truth) + 1]] <<- data.frame(group = grp, protein_id = acc[idx],
                                            gene = genes[idx], true_log2fc = eff,
                                            set = "", stringsAsFactors = FALSE)
}
i_etc <- match(etc, genes); i_tran <- match(tran, genes)
spike(i_etc, "TreatmentA", runif(length(i_etc), 0.8, 1.6))
truth[[length(truth)]]$set <- "electron transport chain (GO:0022900)"
spike(i_tran, "TreatmentB", -runif(length(i_tran), 0.8, 1.6))
truth[[length(truth)]]$set <- "cytoplasmic translation (GO:0002181)"
other <- setdiff(seq_len(n), c(i_etc, i_tran))
for (g in c("TreatmentA", "TreatmentB")) {
  idx <- sample(other, 40)
  spike(idx, g, sample(c(-1, 1), 40, TRUE) * runif(40, 0.8, 1.8))
  truth[[length(truth)]]$set <- "random"
}
truth <- do.call(rbind, truth)

# missingness: detection limit plus a little at random
p_miss <- plogis(-(x - 17.5) * 1.6) + 0.015
x[matrix(runif(length(x)), n) < p_miss] <- NA
lin <- round(2^x, 1)

# ---- example 1: DIA-NN protein-group matrix + sample sheet ----------------
d1 <- file.path(out, "example_diann"); dir.create(d1, recursive = TRUE, showWarnings = FALSE)
pg <- data.frame(Protein.Group = acc, Protein.Ids = acc,
                 Protein.Names = paste0(genes, "_HUMAN"), Genes = genes,
                 First.Protein.Description = unname(desc), check.names = FALSE)
m <- lin; colnames(m) <- paste0("D:\\timsTOF\\2026\\", ss$sample_id, ".d")
write.table(cbind(pg, m), file.path(d1, "report.pg_matrix.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
write.csv(ss, file.path(d1, "samples.csv"), row.names = FALSE)
write.csv(truth, file.path(d1, "truth.csv"), row.names = FALSE)
writeLines(c(
  "study:",
  "  title: \"Example: simulated 3-group cell-line study\"",
  "  description: \"Simulated HEK cell proteomics, two treatments vs control, 5 replicates each (DIA-NN output).\"",
  "  species: human",
  "  tissue: cell line",
  "input:",
  "  matrix: report.pg_matrix.tsv",
  "  format: diann",
  "  sample_sheet: samples.csv",
  "  sample_column: sample_id",
  "design:",
  "  group_column: group",
  "  control: Control",
  "  comparisons: all_vs_control",
  "  covariates: [batch]",
  "preprocessing:",
  "  normalization: median",
  "  imputation: {method: mixed}",
  "statistics:",
  "  method: limma",
  "ai:",
  "  enabled: auto",
  "output:",
  "  dir: results"), file.path(d1, "config.yml"))

# ---- example 2: annotated workbook (core-facility layout) ----------------
# Two groups, sample annotations above the data, and each sample exported as
# a linear quantity plus its log2 (repeated sample names), as some exports do.
if (requireNamespace("openxlsx", quietly = TRUE)) {
  d2 <- file.path(out, "example_annotated"); dir.create(d2, showWarnings = FALSE)
  keep <- ss$group %in% c("Control", "TreatmentA")
  s2 <- ss[keep, ]; l2 <- lin[, keep]
  ids <- rep(s2$sample_id, each = 2)
  vals <- matrix(NA_real_, n, 2 * ncol(l2))
  vals[, seq(1, ncol(vals), 2)] <- l2
  vals[, seq(2, ncol(vals), 2)] <- round(log2(l2), 4)
  vals_chr <- ifelse(is.na(vals), "", format(vals, scientific = FALSE, trim = TRUE))
  top <- rbind(c("Sample", ids),
               c("Condition", rep(s2$group, each = 2)),
               c("Sex", rep(s2$sex, each = 2)),
               c("Batch", rep(s2$batch, each = 2)))
  wb <- rbind(top, cbind(genes, vals_chr))
  openxlsx::write.xlsx(as.data.frame(wb), file.path(d2, "annotated_workbook.xlsx"),
                       colNames = FALSE)
  writeLines(c(
    "study:",
    "  title: \"Example: annotated workbook, 2 groups\"",
    "  species: human",
    "input:",
    "  matrix: annotated_workbook.xlsx",
    "  format: annotated_xlsx",
    "  annotation_rows: [Condition, Sex, Batch]",
    "  duplicate_columns: log2",
    "  is_log2: true",
    "design:",
    "  group_column: Condition",
    "  control: Control",
    "statistics:",
    "  method: limma",
    "ai:",
    "  enabled: false",
    "output:",
    "  dir: results"), file.path(d2, "config.yml"))
  write.csv(truth[truth$group == "TreatmentA", ], file.path(d2, "truth.csv"), row.names = FALSE)
}
cat("Proteins:", n, "| ETC genes:", length(etc), "| translation genes:", length(tran), "\n")
