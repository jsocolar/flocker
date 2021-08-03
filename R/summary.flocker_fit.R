#' Summary method for a flocker_fit object.
#' @param flocker_fit A flocker_fit object produced by \code{flocker(..., visit_constant = F)}.
#' @param ... Other potential arguments.
#' @method summary flocker_fit
#' @export

summary.flocker_fit <- function(flocker_fit) {
  out <- flocker_fit$summary
  class(out) <- "flocker_summary"
}

#' Print a summary for a fitted model represented by a \code{flocker_fit} object.
#' @aliases print.flocker_summary
#' @param flocker_fit An object of class \code{flocker_fit}.
#' @param ... Other potential arguments passed to method \code{summary} of \code{flocker_fit}.
#' @seealso \code{\link{summary.flocker_fit}}
#' @export
print.flocker_fit <- function(flocker_fit, ...) {
  writeLines(paste(summary(flocker_fit, ...), collapse = "\n"))
}