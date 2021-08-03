#' Summary method for a flocker_fit object
#' @param flocker_fit A flocker_fit object produced by \code{flocker(..., visit_constant = F)}
#' @param ... Other potential arguments
#' @method summary flocker_fit
#' @export

summary.flocker_fit <- function(flocker_fit) {
  writeLines(paste(flocker_fit$summary, collapse = "\n"))
}