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
  if (!valid_counts || sum(counts) != nrow(x$data) ||
      !identical(x$data$ff_n_event, rep(as.integer(counts), counts))) {
    stop("Invalid single-season event blocks; reformat the original data.")
  }
  invisible(NULL)
}
