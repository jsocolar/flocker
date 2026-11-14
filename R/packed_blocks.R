# Starts and padding are shared by each level of the nested packing layout.
block_starts <- function(counts) {
  as.integer(cumsum(c(1L, utils::head(counts, -1L))))
}

pad_packed <- function(x, n_row) {
  c(as.integer(x), rep(-99L, n_row - length(x)))
}

# Pack visits in unit_order, retaining mappings to the original input order.
# Inputs are already validated; missing visits must trail within each unit.
# Empty units are opt-in and occupy one placeholder row. Replacement values
# are supplied by the caller independently of the packing order.
pack_unit_block <- function(obs, unit_covs = NULL, event_covs = NULL,
                            allow_empty = FALSE, placeholder_values = NULL,
                            unit_order = seq_len(nrow(obs))) {
  n_unit <- nrow(obs)
  n_event <- as.integer(rowSums(!is.na(obs)))[unit_order]
  if (!allow_empty && any(n_event == 0L)) {
    stop("Empty units are not allowed in this model.")
  }
  n_row_unit <- pmax(1L, n_event)
  unit_start <- block_starts(n_row_unit)
  row_unit <- rep(seq_len(n_unit), n_row_unit)
  row_visit <- sequence(n_row_unit)
  packed_to_original <- unit_order[row_unit] + (row_visit - 1L) * n_unit
  placeholder <- n_event[row_unit] == 0L
  y <- as.integer(obs[packed_to_original])
  y[placeholder] <- -99L
  data <- data.frame(ff_y = y)
  if (!is.null(unit_covs)) {
    data <- cbind(data, unit_covs[unit_order[row_unit], , drop = FALSE])
  }
  for (name in names(event_covs)) {
    values <- event_covs[[name]][packed_to_original]
    missing <- placeholder & is.na(values)
    if (any(missing)) {
      replacement <- placeholder_values[[name]]
      if (length(replacement) != 1L || is.na(replacement)) {
        stop("A replacement value is required for missing placeholder covariate: ", name)
      }
      values[missing] <- replacement
    }
    data[[name]] <- values
  }
  rownames(data) <- NULL
  data$ff_n_unit <- pad_packed(n_unit, nrow(data))
  data$ff_unit_start <- pad_packed(unit_start, nrow(data))
  data$ff_Q <- as.integer(rowSums(obs, na.rm = TRUE) > 0)[unit_order[row_unit]]
  original_to_packed <- matrix(NA_integer_, nrow(obs), ncol(obs), dimnames = dimnames(obs))
  original_to_packed[packed_to_original[!placeholder]] <- which(!placeholder)
  list(data = data, unit_start = unit_start, n_event = n_event, row_unit = row_unit,
       packed_to_original = packed_to_original, placeholder = placeholder,
       original_to_packed = original_to_packed)
}

# Choose detection covariate values for ghost rows representing unvisited units
# in multiseason models. We choose the mean (or for factors the reference level)
# over the entire dataset, independently of series. Only observed visits contribute
# to numeric means.
placeholder_covariate_values <- function(obs, event_covs) {
  lapply(event_covs, function(x) {
    if (is.character(x)) x <- factor(x)
    if (is.factor(x)) return(levels(x)[1L])
    if (is.numeric(x) || is.logical(x)) return(mean(x[!is.na(obs)]))
    stop("Detection covariates must be numeric, logical, character, or factors.")
  })
}

# Use pack_unit_block to pack all units once, then replace the single directory 
# with directories local to each outer block (group or series). Inputs contain 
# retained units only. Input order within each group or series is preserved.
pack_nested_blocks <- function(obs, unit_covs, event_covs, unit_group,
                               outer = c("group", "series"), allow_empty = FALSE,
                               placeholder_values = NULL) {
  outer <- match.arg(outer)
  # Nested models use global factor levels for character detection covariates.
  event_covs <- lapply(event_covs, function(x) if (is.character(x)) factor(x) else x)
  unit_order <- order(unit_group, seq_len(nrow(obs)))
  packed <- pack_unit_block(obs, unit_covs, event_covs, allow_empty,
                             placeholder_values, unit_order)
  group_unit_start <- which(!duplicated(unit_group[unit_order]))
  group_unit_end <- c(group_unit_start[-1L] - 1L, length(unit_order))
  group_start <- packed$unit_start[group_unit_start]
  data <- packed$data
  data$ff_n_unit <- data$ff_unit_start <- rep(-99L, nrow(data))
  for (g in seq_along(group_start)) {
    units <- seq.int(group_unit_start[g], group_unit_end[g])
    directory <- group_start[g] + seq_along(units) - 1L
    data$ff_n_unit[group_start[g]] <- length(units)
    data$ff_unit_start[directory] <- packed$unit_start[units] - group_start[g] + 1L
  }
  data[[paste0("ff_n_", outer)]] <- pad_packed(length(group_start), nrow(data))
  data[[paste0("ff_", outer, "_start")]] <- pad_packed(group_start, nrow(data))
  # Preserve the public column order: outer directory precedes unit directory.
  directory <- c("ff_n_unit", "ff_unit_start", "ff_Q")
  packed$data <- data[c(setdiff(names(data), directory), directory)]
  packed$unit_order <- unit_order
  packed$group_start <- group_start
  packed$unit_to_packed <- integer(nrow(obs))
  packed$unit_to_packed[unit_order] <- packed$unit_start
  names(packed$unit_to_packed) <- rownames(obs)
  packed
}
