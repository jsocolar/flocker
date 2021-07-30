summary.flocker_fit <- function(flocker_fit) {
  writeLines(paste(flocker_fit$summary, collapse = "\n"))
}