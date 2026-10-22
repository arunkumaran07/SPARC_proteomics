#' Example datasets shipped with the package
#'
#' Copies an example dataset to a writable folder and returns the path to its
#' configuration, ready for [run_pipeline()]. The data are **simulated**, but
#' use real human gene symbols and UniProt accessions so that annotation and
#' Reactome enrichment behave as on real data.
#'
#' \describe{
#'   \item{`"diann"`}{DIA-NN `report.pg_matrix.tsv` (1,800 proteins), a sample
#'     sheet, and a config. Three groups (Control, TreatmentA, TreatmentB) x 5
#'     replicates, with sex and batch columns. TreatmentA raises proteins
#'     annotated to the electron transport chain (GO:0022900); TreatmentB
#'     lowers proteins annotated to cytoplasmic translation (GO:0002181); each
#'     treatment also changes 40 random proteins.}
#'   \item{`"annotated"`}{The same proteins as a core-facility style Excel
#'     workbook: sample annotations above the data, Control vs TreatmentA, and
#'     each sample exported twice (linear and log2) under the same name.}
#' }
#' Both folders include `truth.csv`, listing every planted change.
#'
#' @param name `"diann"` or `"annotated"`.
#' @param dir Folder to copy the example into (created if needed).
#' @return Path to the example's `config.yml`.
#' @export
#' @examples
#' cfg <- example_data("diann", tempfile("example"))
#' list.files(dirname(cfg))
#' \donttest{
#' res <- run_pipeline(cfg)
#' res$de$summary
#' }
example_data <- function(name = c("diann", "annotated"), dir = "sparc_example") {
  name <- match.arg(name)
  src <- system.file("extdata", paste0("example_", name), package = "sparcproteomics")
  if (!nzchar(src)) stop("Example data not found; reinstall the package.", call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(src, full.names = TRUE), dir, overwrite = TRUE)
  if (!all(ok)) stop("Could not copy the example into ", dir, call. = FALSE)
  file.path(dir, "config.yml")
}
