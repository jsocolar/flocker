#' Summary method for a flocker_fit object
#' @param flocker_fit A flocker_fit object produced by \code{flocker(..., visit_constant = F)}
#' @return Print a summary output for the fit.
#' @export

summary.flocker_fit <- function(flocker_fit) {
  writeLines(paste(flocker_fit$summary, collapse = "\n"))
}