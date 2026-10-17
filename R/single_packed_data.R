# Validate the new rep-varying fitting contract before passing rows to brms.
validate_single_packed_data <- function(x) {
  metadata <- get_flocker_metadata(x)
  if (!identical(metadata$schema_version, 2L)) {
    stop("Incompatible single-season packing schema. Pre-2.0 objects require a compatible 1.x installation or reformatting and refitting.")
  }
  if (!identical(metadata$data_role, "model_data")) {
    stop("Fitting requires model_data; new_data cannot be used for fitting.")
  }
  counts <- metadata$n_event
  valid_counts <- is.numeric(counts) && length(counts) > 0L &&
    all(is.finite(counts)) && all(counts >= 1 & counts == floor(counts))
  if (!valid_counts || sum(counts) != nrow(x$data)) {
    stop("Invalid single-season event blocks; reformat the original data.")
  }
  n_unit <- length(counts)
  starts <- cumsum(c(1L, head(as.integer(counts), -1L)))
  if (!identical(metadata$unit_start, starts) ||
      !identical(x$data$ff_n_unit, c(n_unit, rep(-99L, nrow(x$data) - 1L))) ||
      !identical(x$data$ff_unit_start, c(starts, rep(-99L, nrow(x$data) - n_unit)))) {
    stop("Invalid single-season event blocks; reformat the original data.")
  }
  invisible(NULL)
}
