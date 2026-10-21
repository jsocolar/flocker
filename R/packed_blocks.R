# Starts and padding are shared by each level of the nested packing layout.
block_starts <- function(counts) {
  as.integer(cumsum(c(1L, head(counts, -1L))))
}

pad_packed <- function(x, n_row) {
  c(as.integer(x), rep(-99L, n_row - length(x)))
}

# Pack visits within units and units within outer groups. Inputs have one row
# per retained unit, in the caller's desired within-group order. Model-specific
# validation and selection of retained units belong to the caller. In particular,
# multiseason models allow empty units; two-level models need not allow them.
pack_nested_blocks <- function(obs, unit_covs, event_covs, unit_group) {
  unit_order <- order(unit_group, seq_len(nrow(obs)))
  group <- unit_group[unit_order]
  n_event <- as.integer(rowSums(!is.na(obs)))[unit_order]
  n_row_unit <- pmax(1L, n_event)
  unit_start <- block_starts(n_row_unit)
  group_unit_start <- which(!duplicated(group))
  group_start <- unit_start[group_unit_start]
  row_unit <- rep(seq_along(unit_order), n_row_unit)
  row_visit <- sequence(n_row_unit)
  packed_to_original <- unit_order[row_unit] + (row_visit - 1L) * nrow(obs)
  placeholder <- n_event[row_unit] == 0L
  y <- as.integer(obs[packed_to_original])
  y[placeholder] <- -99L
  data <- data.frame(ff_y = y)
  if (!is.null(unit_covs)) {
    data <- cbind(data, unit_covs[unit_order[row_unit], , drop = FALSE])
  }
  for (name in names(event_covs)) {
    x <- event_covs[[name]]
    # Character arrays use the same default level ordering as factor().
    if (is.character(x)) x <- factor(x)
    values <- x[packed_to_original]
    missing <- placeholder & is.na(values)
    if (any(missing)) {
      if (is.factor(x)) {
        values[missing] <- levels(x)[1L]
      } else if (is.numeric(x) || is.logical(x)) {
        values[missing] <- mean(x[!is.na(obs)])
      } else {
        stop("Detection covariates must be numeric, logical, character, or factors.")
      }
    }
    data[[name]] <- values
  }
  rownames(data) <- NULL
  original_to_packed <- matrix(NA_integer_, nrow(obs), ncol(obs), dimnames = dimnames(obs))
  original_to_packed[packed_to_original[!placeholder]] <- which(!placeholder)
  # Each group carries its own unit directory in its leading rows. Unit
  # starts are relative to the group, so an extracted block needs no rebasing.
  local_n_unit <- local_unit_start <- local_Q <- rep(-99L, nrow(data))
  Q <- as.integer(rowSums(obs, na.rm = TRUE)[unit_order] > 0)
  group_unit_end <- c(group_unit_start[-1L] - 1L, length(unit_start))
  for (g in seq_along(group_start)) {
    units <- seq.int(group_unit_start[g], group_unit_end[g])
    directory <- group_start[g] + seq_along(units) - 1L
    local_n_unit[group_start[g]] <- length(units)
    local_unit_start[directory] <- unit_start[units] - group_start[g] + 1L
    local_Q[directory] <- Q[units]
  }
  list(data = data, unit_order = unit_order, group_start = as.integer(group_start),
       unit_start = unit_start, n_event = n_event, row_unit = row_unit,
       packed_to_original = packed_to_original, placeholder = placeholder,
       original_to_packed = original_to_packed,
       n_unit = local_n_unit, local_unit_start = local_unit_start, Q = local_Q)
}
