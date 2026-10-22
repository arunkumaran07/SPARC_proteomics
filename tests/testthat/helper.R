withr_tempdir <- function() {
  d <- tempfile("sp_test")
  dir.create(d)
  d
}

# Fraction of true changes recovered and false discovery proportion
recovery <- function(res, truth_path) {
  truth <- utils::read.csv(truth_path)
  out <- lapply(names(res$de$tables), function(nm) {
    g <- sub("_vs_Control$", "", nm)
    tp <- truth$protein_id[truth$group == g]
    sig <- res$de$tables[[nm]]$protein_id[res$de$tables[[nm]]$significant]
    c(power = mean(tp %in% sig), fdp = if (length(sig)) mean(!sig %in% tp) else 0)
  })
  do.call(rbind, out)
}
