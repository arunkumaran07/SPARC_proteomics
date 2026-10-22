.pkg_state <- new.env(parent = emptyenv())
.pkg_state$log <- character(0)
.pkg_state$verbose <- TRUE

# Timestamped progress message, also kept for the run log written by
# run_pipeline().
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "[%H:%M:%S] "), paste0(..., collapse = ""))
  .pkg_state$log <- c(.pkg_state$log, msg)
  if (isTRUE(.pkg_state$verbose)) message(msg)
  invisible(msg)
}

log_reset <- function() .pkg_state$log <- character(0)

#' Messages logged during the current session's analysis
#'
#' @return A character vector of timestamped messages.
#' @export
get_run_log <- function() .pkg_state$log

.need <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is needed ", why, ". Install it with ",
         if (pkg %in% c("impute", "clusterProfiler", "ReactomePA") ||
             grepl("^org\\.", pkg))
           paste0("BiocManager::install(\"", pkg, "\")")
         else paste0("install.packages(\"", pkg, "\")"),
         ".", call. = FALSE)
}
