#' Summary method for a flocker_fit object.
#' @param object A flocker_fit object produced by \code{flocker(..., visit_constant = F)}.
#' @param ... Other potential arguments.
#' @method summary flocker_fit
#' @export

summary.flocker_fit <- function(object, ...) {
  out <- object$summary
  class(out) <- "flocker_summary"
  out
}

#' Print a summary for a fitted model represented by a \code{flocker_fit} object.
#' @aliases print.flocker_summary
#' @param x An object of class \code{flocker_fit}.
#' @param ... Other potential arguments passed to method \code{summary} of \code{flocker_fit}.
#' @seealso \code{\link{summary.flocker_fit}}
#' @export
print.flocker_fit <- function(x, ...) {
  writeLines(paste(summary(x, ...), collapse = "\n"))
}
