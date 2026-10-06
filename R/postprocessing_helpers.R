#' Prepare probability components for two-level post-processing
#' @noRd
prepare_twolevel_postprocessing <- function(
    flocker_fit, lik_type, draw_ids, new_data, allow_new_levels,
    sample_new_levels, include_detection
    ) {
  data_object <- if(is.null(new_data)) flocker_fit else new_data
  the_data <- data_object$data
  metadata <- get_flocker_metadata(data_object)

  unit_lps <- fitted_flocker(
    flocker_fit, components = c("occ", "Omega"), draw_ids = draw_ids,
    new_data = new_data, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels, response = TRUE,
    unit_level = TRUE
  )

  if(include_detection) {
    gp <- get_positions(data_object)
    obs <- new_array(gp, the_data$ff_y[gp])
    det_lps <- fitted_flocker(
      flocker_fit, components = "det", draw_ids = draw_ids,
      new_data = new_data, allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels, response = TRUE,
      unit_level = FALSE
    )
  } else {
    obs <- NULL
    det_lps <- NULL
  }

  format_twolevel_postprocessing(
    unit_lps, det_lps, obs, lik_type, the_data, metadata
  )
}

#' Put two-level predictions in a common unit-by-draw representation
#' @noRd
format_twolevel_postprocessing <- function(
    unit_lps, det_lps, obs, lik_type, flocker_data_data, flocker_metadata
    ) {
  n_unit <- flocker_data_data$ff_n_unit[1]
  unit_rows <- seq_len(n_unit)
  Omega <- group_level_Omega(
    unit_lps, lik_type, flocker_data_data, flocker_metadata
  )
  n_group <- nrow(Omega)
  n_draw <- ncol(Omega)
  group_known_present <- flocker_data_data$ff_group_known_present[
    seq_len(n_group)
  ]

  if(lik_type == "twolevel_single") {
    psi_all <- matrix(unit_lps$linpred_occ, nrow = n_unit)
    theta_all <- if(is.null(det_lps)) NULL else det_lps$linpred_det
    obs_use <- obs
    packed_group_id <- get_unit_group(flocker_data_data)
    group_id <- integer(n_unit)
    group_id[flocker_metadata$unit_order] <- packed_group_id
    unit_site <- NULL
    n_site <- NULL
  } else {
    assertthat::assert_that(lik_type == "augmented")
    unit_site <- flocker_metadata$unit_site
    group_id <- get_unit_group(flocker_data_data)
    n_site <- max(unit_site)
    psi_array <- array(
      unit_lps$linpred_occ, dim = c(n_site, n_group, n_draw)
    )
    psi_all <- matrix(NA_real_, nrow = n_unit, ncol = n_draw)
    for(i in unit_rows) {
      psi_all[i, ] <- psi_array[unit_site[i], group_id[i], ]
    }

    if(is.null(det_lps)) {
      theta_all <- NULL
      obs_use <- NULL
    } else {
      n_visit <- dim(obs)[2]
      theta_array <- array(
        det_lps$linpred_det, dim = c(n_site, n_visit, n_group, n_draw)
      )
      theta_all <- array(NA_real_, dim = c(n_unit, n_visit, n_draw))
      obs_use <- matrix(NA_real_, nrow = n_unit, ncol = n_visit)
      for(i in unit_rows) {
        theta_all[i, , ] <- theta_array[unit_site[i], , group_id[i], ]
        obs_use[i, ] <- obs[unit_site[i], , group_id[i]]
      }
    }
  }

  list(
    psi = psi_all,
    theta = theta_all,
    Omega = Omega,
    group_id = group_id,
    group_known_present = group_known_present,
    obs = obs_use,
    group_names = get_level2_group_names(
      flocker_data_data, flocker_metadata, lik_type
    ),
    unit_site = unit_site,
    n_site = n_site
  )
}

#' Compute occupancy-model emission likelihoods by unit and draw
#' @noRd
occupancy_emission_components <- function(psi, theta, obs) {
  n_unit <- nrow(psi)
  n_draw <- ncol(psi)
  el_0 <- el_1 <- matrix(NA_real_, nrow = n_unit, ncol = n_draw)
  for(i in seq_len(n_draw)) {
    el_0[, i] <- emission_likelihood(0, obs, theta[, , i])
    el_1[, i] <- emission_likelihood(1, obs, theta[, , i])
  }
  list(
    unavailable = el_0,
    available = el_1,
    unit_lik_available = (1 - psi) * el_0 + psi * el_1
  )
}

#' Get level-two group names in internal group order
#' @noRd
get_level2_group_names <- function(flocker_data_data, flocker_metadata,
                                   data_type) {
  group_names <- flocker_metadata$level2_group_names
  n_group <- flocker_data_data$ff_n_group[1]

  if (is.null(group_names)) {
    if (data_type == "twolevel_single") {
      group_names <- as.character(
        flocker_data_data[[flocker_metadata$level2_group]][seq_len(n_group)]
      )
    } else {
      assertthat::assert_that(data_type == "augmented")
      group_names <- paste0("species", seq_len(n_group))
    }
  }

  assertthat::assert_that(length(group_names) == n_group)
  as.character(group_names)
}

#' Recover one Omega value per level-two group from fitted_flocker output
#' @noRd
group_level_Omega <- function(fitted_output, data_type, flocker_data_data,
                              flocker_metadata) {
  assertthat::assert_that("linpred_Omega" %in% names(fitted_output))
  unit_level <- attr(fitted_output, "unit_level")
  assertthat::assert_that(is_one_logical(unit_level))
  Omega <- fitted_output$linpred_Omega
  n_group <- flocker_data_data$ff_n_group[1]
  if(data_type == "twolevel_single") {
    n_unit <- flocker_data_data$ff_n_unit[1]
    if(unit_level) {
      unit_Omega <- matrix(Omega, nrow = n_unit)
    } else if(length(dim(Omega)) == 2) {
      unit_Omega <- matrix(Omega[, 1], ncol = 1)
    } else {
      unit_Omega <- matrix(Omega[, 1, , drop = FALSE], nrow = n_unit)
    }
    group_representatives <- flocker_metadata$unit_order[seq_len(n_group)]
    out <- unit_Omega[group_representatives, , drop = FALSE]
  } else {
    assertthat::assert_that(data_type == "augmented")
    if(unit_level) {
      if(length(dim(Omega)) == 2) {
        out <- matrix(Omega[1, ], ncol = 1)
      } else {
        out <- Omega[1, , , drop = FALSE]
      }
    } else {
      if(length(dim(Omega)) == 3) {
        out <- matrix(Omega[1, 1, ], ncol = 1)
      } else {
        out <- Omega[1, 1, , , drop = FALSE]
      }
    }
  }
  matrix(out, nrow = n_group)
}
