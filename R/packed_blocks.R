# Starts and padding are shared by each level of the nested packing layout.
block_starts <- function(counts) {
  as.integer(cumsum(c(1L, head(counts, -1L))))
}

pad_packed <- function(x, n_row) {
  c(as.integer(x), rep(-99L, n_row - length(x)))
}

# Pack visits within one ordered block of units. Public formatters validate
# outcomes, covariates, and trailing missing visits before calling this helper.
# Empty units are opt-in and occupy one placeholder row. Replacement values
# are supplied by the caller so splitting into blocks cannot change imputation.
pack_unit_block <- function(obs, unit_covs = NULL, event_covs = NULL,
                            allow_empty = FALSE, placeholder_values = NULL) {
  n_unit <- nrow(obs)
  n_event <- as.integer(rowSums(!is.na(obs)))
  if (!allow_empty && any(n_event == 0L)) {
    stop("Empty units are not allowed in this model.")
  }
  n_row_unit <- pmax(1L, n_event)
  unit_start <- block_starts(n_row_unit)
  row_unit <- rep(seq_len(n_unit), n_row_unit)
  row_visit <- sequence(n_row_unit)
  packed_to_original <- row_unit + (row_visit - 1L) * n_unit
  placeholder <- n_event[row_unit] == 0L
  y <- as.integer(obs[packed_to_original])
  y[placeholder] <- -99L
  data <- data.frame(ff_y = y)
  if (!is.null(unit_covs)) {
    data <- cbind(data, unit_covs[row_unit, , drop = FALSE])
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
  data$ff_Q <- as.integer(rowSums(obs, na.rm = TRUE) > 0)[row_unit]
  original_to_packed <- matrix(NA_integer_, nrow(obs), ncol(obs), dimnames = dimnames(obs))
  original_to_packed[packed_to_original[!placeholder]] <- which(!placeholder)
  list(data = data, unit_start = unit_start, n_event = n_event, row_unit = row_unit,
       packed_to_original = packed_to_original, placeholder = placeholder,
       original_to_packed = original_to_packed)
}

# Choose replacements over the entire dataset, before splitting into series.
# Only observed visits contribute to numeric means. Factor levels retain their
# global ordering, even if a level is absent from an individual series.
placeholder_covariate_values <- function(obs, event_covs) {
  lapply(event_covs, function(x) {
    if (is.character(x)) x <- factor(x)
    if (is.factor(x)) return(levels(x)[1L])
    if (is.numeric(x) || is.logical(x)) return(mean(x[!is.na(obs)]))
    stop("Detection covariates must be numeric, logical, character, or factors.")
  })
}

# Pack each outer group using the same one-level layout as ordinary single-
# season models. Inputs contain retained units only, in within-group order.
# outer names the directory columns; allow_empty controls placeholder units.
pack_nested_blocks <- function(obs, unit_covs, event_covs, unit_group,
                               outer = c("group", "series"), allow_empty = FALSE,
                               placeholder_values = NULL) {
  outer <- match.arg(outer)
  # Normalize character levels before splitting, including levels absent from
  # individual groups. Single-season data retain their existing column types.
  event_covs <- lapply(event_covs, function(x) {
    if (is.character(x)) {
      values <- factor(x)
      dim(values) <- dim(x)
      return(values)
    }
    x
  })
  unit_order <- order(unit_group, seq_len(nrow(obs)))
  group <- unit_group[unit_order]
  group_unit_start <- which(!duplicated(group))
  group_unit_end <- c(group_unit_start[-1L] - 1L, length(unit_order))
  blocks <- lapply(seq_along(group_unit_start), function(g) {
    units <- unit_order[seq.int(group_unit_start[g], group_unit_end[g])]
    indices <- base::outer(units, (seq_len(ncol(obs)) - 1L) * nrow(obs), `+`)
    events <- lapply(event_covs, function(x) {
      # Linear indexing preserves factor classes and their global levels.
      values <- x[as.vector(indices)]
      dim(values) <- dim(indices)
      values
    })
    pack_unit_block(obs[units, , drop = FALSE],
                    if (is.null(unit_covs)) NULL else unit_covs[units, , drop = FALSE],
                    events, allow_empty, placeholder_values)
  })
  group_start <- block_starts(vapply(blocks, function(x) nrow(x$data), integer(1)))
  data <- do.call(rbind, lapply(blocks, `[[`, "data"))
  rownames(data) <- NULL
  data[[paste0("ff_n_", outer)]] <- pad_packed(length(blocks), nrow(data))
  data[[paste0("ff_", outer, "_start")]] <- pad_packed(group_start, nrow(data))
  # Preserve the public column order: outer directory precedes unit directory.
  directory <- c("ff_n_unit", "ff_unit_start", "ff_Q")
  data <- data[c(setdiff(names(data), directory), directory)]
  unit_start <- as.integer(unlist(lapply(seq_along(blocks), function(g) {
    blocks[[g]]$unit_start + group_start[g] - 1L
  })))
  n_event <- unlist(lapply(blocks, `[[`, "n_event"))
  row_unit <- rep(seq_along(unit_order), pmax(1L, n_event))
  packed_to_original <- unit_order[row_unit] +
    (sequence(pmax(1L, n_event)) - 1L) * nrow(obs)
  placeholder <- unlist(lapply(blocks, `[[`, "placeholder"))
  original_to_packed <- matrix(NA_integer_, nrow(obs), ncol(obs), dimnames = dimnames(obs))
  original_to_packed[packed_to_original[!placeholder]] <- which(!placeholder)
  unit_to_packed <- integer(nrow(obs))
  unit_to_packed[unit_order] <- unit_start
  names(unit_to_packed) <- rownames(obs)
  list(data = data, unit_order = unit_order, group_start = group_start,
       unit_start = unit_start, n_event = n_event, row_unit = row_unit,
       packed_to_original = packed_to_original, placeholder = placeholder,
       original_to_packed = original_to_packed, unit_to_packed = unit_to_packed)
}
