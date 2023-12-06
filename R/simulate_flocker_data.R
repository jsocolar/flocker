#' Simulate data for use with \code{make_flocker_data()} and downstream 
#' functions. 
#' 
#' Data will be simulated with one unit covariate that affects all 
#' relevant terms, one event covariate that affects detection (unless 
#' `rep_constant` is `TRUE`), and one grouping factor representing species
#' with correlated effects on all terms.
#' @param n_rep number of replicate visits to simulate per closure unit
#' @param n_pt number of points to simulate. The number of units for single-
#'   season models will be `n_pt*n_sp`. The number of units for multi-season
#'   models will be `n_pt*n_sp*n_season`.
#' @param n_sp number of levels to include in random effect. For compatibility
#'   with multispecies models where the random effect represents species,
#'   the data get expanded such that there's a row (closure-unit) for each
#'   combination of sampling point and effect level (i.e. species)
#' @param n_season Number of seasons desired. 1 yields data for a single-season 
#'   model; all other positive integers yield data for multiseason models.
#' @param multiseason if n_season is NULL, must be NULL. Otherwise, one of
#'   "colex" or "autologistic".
#' @param multi_init if n_season is NULL, must be NULL. Otherwise, one of 
#'   "explicit" or "equilibrium".
#' @param augmented logical. If `TRUE` data will be formatted for an augmented
#'   model, which requires that `n_season == 1`. All never-observed 
#'   species will be trimmed out of the data, and the default parameters will be
#'   modified to increase random effect variances for the detection and occupancy,
#'   intercepts and to decrease random effect variances for detection slopes, thus
#'   encouraging the existence of never-observed species. Furthermore, the data will
#'   be simulated without any covariate influence on occupancy.
#' @param rep_constant logical: create data with unit covariates only (TRUE) 
#'   or data that includes event covariates (FALSE)
#' @param params a named list containing of parameter values to use in simulation.
#'   Any required parameters whose names are not in this list will be assigned their
#'   default values. To see the parameter names and structures required, run with 
#'   `params = NULL` (the default) and examine the `$params` element of the output.
#' @param covariates a dataframe of covariate values to use in simulation, or
#'   NULL to simulate values. To see the covariate names and structures required,
#'   run with `covariates = NULL` (the default) and examine the `$covariates` element
#'   of the output.
#' @param seed random seed. NULL uses (and updates) the existing RNG state. Other values
#'   do not update the global RNG state.
#' @param ragged_rep logical: create data with variable (TRUE) or constant 
#'   (FALSE) numbers of visits per unit.  If TRUE, approximately half of units 
#'   will be missing approximately half of `n_rep` visits. Intended primarily for 
#'   development purposes (bug-checking models with missing data).
#' @param missing_seasons logical; relevant only if n_season is greater than 1. 
#'   create data with variable (TRUE) or constant (FALSE) numbers of seasons per 
#'   series (TRUE). If TRUE, approximately half of series will be missing their
#'   even-numbered seasons.
#' @return A named list with the observation matrix/array ($obs), the unit covariate 
#'   dataframe(s) ($unit_covs), the event covariate list ($event_covs), the parameters
#'   used in simulation ($params) and the covariate list used in simulation ($covariates). 
#'   If rep_constant is TRUE, then $event_covs will be NULL. 
#' @export
#' @examples 
#' simulate_flocker_data()
simulate_flocker_data <- function(
  n_rep = 4,
  n_pt = 50,
  n_sp = 30,
  n_season = 1,
  multiseason = NULL,
  multi_init = NULL,
  augmented = FALSE,
  rep_constant = FALSE,
  params = NULL,
  covariates = NULL,
  seed = 123,
  ragged_rep = FALSE,
  missing_seasons = FALSE) {
  
  # Input checking
  msg_str <- " is not a single positive integer"
  assertthat::assert_that(is_one_pos_int(n_rep, 1), msg = paste0("n_rep", msg_str, " greater than one"))
  assertthat::assert_that(is_one_pos_int(n_pt, 0), msg = paste0("n_pt", msg_str))
  assertthat::assert_that(is_one_pos_int(n_sp, 0), msg = paste0("n_sp", msg_str))
  assertthat::assert_that(is_one_pos_int(n_season, 0), msg = paste0("n_season", msg_str))
  assertthat::assert_that(
    is.null(multiseason) | isTRUE(multiseason %in% c("colex", "autologistic")),
    msg = "multiseason must be NULL, 'colex', or 'autologistic'"
  )
  assertthat::assert_that(
    is.null(multi_init) | isTRUE(multi_init %in% c("explicit", "equilibrium")),
    msg = "multiseason must be NULL, 'explicit', or 'equilibrium'"
  )
  assertthat::assert_that(is_one_logical(augmented), msg = "augmented must be a single logical")
  assertthat::assert_that(is_one_logical(rep_constant), msg = "rep_constant must be a single logical")
  assertthat::assert_that(is_one_logical(ragged_rep), msg = "ragged_rep must be a single logical")
  assertthat::assert_that(is_one_logical(missing_seasons), msg = "missing_seasons must be a single logical")
  if(n_season == 1){
    assertthat::assert_that(
      is.null(multiseason) & is.null(multi_init),
      msg = "multiseason and multi_init must be NULL when n_season is 1"
    )
    assertthat::assert_that(
      !missing_seasons,
      msg = "missing_seasons must be FALSE when n_season is 1"
    )
  } else {
    assertthat::assert_that(
      !is.null(multiseason) & !is.null(multi_init),
      msg = "multiseason and multi_init must be provide when n_season is greater than 1"
    )
    assertthat::assert_that(
      !augmented, 
      msg = "augmented must be FALSE when n_season is greater than 1"
    )
  }
  
  if (is.null(seed)){
    out <- sfd(n_rep, n_pt, n_sp, n_season, multiseason, multi_init, augmented,
               rep_constant, params, covariates, ragged_rep, missing_seasons)
  } else {
      out <- withr::with_seed(
        seed,
        sfd(n_rep, n_pt, n_sp, n_season, multiseason, multi_init, augmented,
            rep_constant, params, covariates, ragged_rep, missing_seasons)
      )
  }
  out
}
  
#' util for creating example data
#' @inherit simulate_flocker_data
#' @noRd
sfd <- function(
    n_rep,
    n_pt,
    n_sp,
    n_season,
    multiseason,
    multi_init,
    augmented,
    rep_constant,
    params,
    covariates,
    ragged_rep,
    missing_seasons
) {
  
  if (is.null(params)) {
    params <- list()
  }
  
  # list out coefficients used by model type
  coef_names <- c("det_intercept", "det_slope_unit")
  if (!rep_constant) {coef_names <- c(coef_names, "det_slope_visit")}
  if (n_season == 1 | isTRUE(multi_init == "explicit")) {
    coef_names <- c(coef_names, c("occ_intercept"))
    if (!augmented) {coef_names <- c(coef_names, "occ_slope_unit")}
  }
  if (n_season > 1) {
    coef_names <- c(coef_names, c("col_intercept", "col_slope_unit"))
  }
  if (isTRUE(multiseason == "colex")) {
    coef_names <- c(coef_names, c("ex_intercept", "ex_slope_unit"))
  } else if (isTRUE(multiseason == "autologistic")){
    coef_names <- c(coef_names, c("auto_intercept", "auto_slope_unit"))
  }
  
  # If coefs not specified directly, simulate from the prior (either passed via
  # params or otherwise the default prior).
  n_coef <- length(coef_names)
  np <- names(params)
  if (!("coefs" %in% np)) { # coef values not specified directly
    if (!("coef_means" %in% np)) {
      params$coef_means <- rep(0, n_coef)
    }
    if (!("Sigma" %in% np)) {
      Sigma <- matrix(.5, nrow = n_coef, ncol = n_coef)
      diag(Sigma) <- 1
      if(augmented){ # fiddle with variances to encourage never-detected species in the data
        Sigma[which(coef_names == "det_intercept"), ] <- 1.7 * Sigma[which(coef_names == "det_intercept"), ]
        Sigma[, which(coef_names == "det_intercept")] <- 1.7 * Sigma[, which(coef_names == "det_intercept")]
        Sigma[which(coef_names == "occ_intercept"), ] <- 1.7 * Sigma[which(coef_names == "occ_intercept"), ]
        Sigma[, which(coef_names == "occ_intercept")] <- 1.7 * Sigma[, which(coef_names == "occ_intercept")]
        Sigma[which(coef_names == "det_slope_unit"), ] <- .3 * Sigma[which(coef_names == "det_slope_unit"), ]
        Sigma[, which(coef_names == "det_slope_unit")] <- .3 * Sigma[, which(coef_names == "det_slope_unit")]
        if(!rep_constant) {
          Sigma[which(coef_names == "det_slope_visit"), ] <- .3 * Sigma[which(coef_names == "det_slope_visit"), ]
          Sigma[, which(coef_names == "det_slope_visit")] <- .3 * Sigma[, which(coef_names == "det_slope_visit")]
        }
      }
      params$Sigma <- Sigma
    }
    if (n_sp == 1) {
      params$coefs <- as.data.frame(
        t(as.matrix(MASS::mvrnorm(n_sp, params$coef_means, params$Sigma)))
      )
    } else {
      params$coefs <- as.data.frame(
        MASS::mvrnorm(n_sp, params$coef_means, params$Sigma)
      )
    }

    names(params$coefs) <- coef_names
  }
  assertthat::assert_that(
    identical(names(params$coefs), coef_names)
  )
  coefs <- params$coefs
  
  # Create covariates (currently assume that unit covariates are the same
  # across all timesteps in multiseason models)
  unit_backbone <- expand.grid(
    species = factor(paste0("sp_", 1:n_sp), levels = paste0("sp_", 1:n_sp)),
    id_point = 1:n_pt
    )
  if(is.null(covariates)) {
    covariates <- list(
      uc1 = stats::rnorm(n_pt)[unit_backbone$id_point]
    )
  }
  
  unit_backbone$uc1 <- covariates$uc1
  
  if (n_season > 1) {
    unit_backbone$col <- coefs$col_intercept[as.integer(unit_backbone$species)] +
      coefs$col_slope_unit[as.integer(unit_backbone$species)] * unit_backbone$uc1
    if (multiseason == "colex") {
      unit_backbone$ex <- coefs$ex_intercept[as.integer(unit_backbone$species)] +
        coefs$ex_slope_unit[as.integer(unit_backbone$species)] * unit_backbone$uc1
    } else if (multiseason == "autologistic") {
      unit_backbone$ex <- -(unit_backbone$col + 
        coefs$auto_intercept[as.integer(unit_backbone$species)] +
        coefs$auto_slope_unit[as.integer(unit_backbone$species)] * unit_backbone$uc1)
    }
  }
  
  if (augmented) {
    unit_backbone$occ <- coefs$occ_intercept[as.integer(unit_backbone$species)]
  } else if(n_season == 1 | isTRUE(multi_init == "explicit")){
    unit_backbone$occ <- coefs$occ_intercept[as.integer(unit_backbone$species)] +
      coefs$occ_slope_unit[as.integer(unit_backbone$species)] * unit_backbone$uc1
  } else if(isTRUE(multi_init == "equilibrium")) {
    unit_backbone$occ <- boot::logit(
      boot::inv.logit(unit_backbone$col)/
        (boot::inv.logit(unit_backbone$col) + boot::inv.logit(unit_backbone$ex))
    )
  }
  
  true_Z <- cbind(
    unit_backbone,
    data.frame(
      id_season = 1,
      true_Z = stats::rbinom(n_sp * n_pt, 1, boot::inv.logit(unit_backbone$occ))
    )
  )
  
  if (n_season > 1) {
    for(i in 2:n_season) {
      psi_i <- true_Z$true_Z[true_Z$id_season == (i - 1)] * boot::inv.logit(-unit_backbone$ex) +
        (!true_Z$true_Z[true_Z$id_season == (i - 1)]) * boot::inv.logit(unit_backbone$col)
      true_Z <- 
        rbind(
          true_Z,
          cbind(
            unit_backbone,
            data.frame(
              id_season = i, 
              true_Z = stats::rbinom(n_pt * n_sp, 1, psi_i)
            )
          )
        )
    }
  }
  
  visit_backbone <- expand.grid(
    species = factor(paste0("sp_", 1:n_sp), levels = paste0("sp_", 1:n_sp)),
    id_point = 1:n_pt,
    id_season = 1:n_season,
    id_rep = 1:n_rep
  )
  
  visit_backbone$rowid <- seq_len(nrow(visit_backbone))
  vf1 <- merge(
    visit_backbone, 
    true_Z, 
    by = c("species", "id_point", "id_season"),
    all = TRUE
  )
  visit_full <- vf1[order(vf1$rowid), ]
  
  visit_full$logit_det <- 
    coefs$det_intercept[as.integer(visit_full$species)] +
    coefs$det_slope_unit[as.integer(visit_full$species)] * visit_full$uc1
  
  if (!rep_constant) {
    if (!("ec1" %in% names(covariates))) {
      covariates$ec1 <- stats::rnorm(n_pt*n_season*n_rep)[ # we don't allow ec1 to vary by species;
        # this allows the covariate to play nicely with augmented models.
        
        # I think that visit_full is already guaranteed to be blocked by species with everything
        # else in consistent order, but this should ensure that the ordering stays correct no 
        # matter what.
        as.integer(factor( 
          paste(
            visit_full$id_point, 
            visit_full$id_season, 
            visit_full$id_rep, 
            sep = "__"
          )
        ))
      ]
    }
    visit_full$ec1 <- covariates$ec1
    visit_full$logit_det <- visit_full$logit_det + 
      coefs$det_slope_visit[as.integer(visit_full$species)] * visit_full$ec1
  }
  
  visit_full$obs <- visit_full$true_Z *
    stats::rbinom(nrow(visit_full), 1, boot::inv.logit(visit_full$logit_det))
  
  if (ragged_rep) {
    visit_full$season_pt <- interaction(visit_full$id_point, visit_full$id_season)
    season_pt_missing <- sample(unique(visit_full$season_pt), floor(length(unique(visit_full$season_pt)) / 2))
    rows_missing <- which((visit_full$season_pt %in% season_pt_missing) & (visit_full$id_rep > n_rep / 2))
    visit_full$obs[rows_missing] <- NA
    if(!rep_constant){
      visit_full$ec1[rows_missing] <- NA
    }
  }
  
  if (missing_seasons) {
    pt_missing <- sample(unique(visit_full$id_point), floor(length(unique(visit_full$id_point)) / 2))
    rows_missing <- which((visit_full$id_point %in% pt_missing) & (visit_full$id_season %% 2 == 0))
    visit_full$obs[rows_missing] <- visit_full$ec1[rows_missing] <- NA
  }
  
  # format data for return
  if (! augmented & n_season == 1) {
    visit_full$id_unit <- interaction(visit_full$id_point, visit_full$species)
    unit_backbone$id_unit <- interaction(unit_backbone$id_point, unit_backbone$species)
    obs <- t(utils::unstack(visit_full[c("obs", "id_unit")], obs ~ id_unit))
    if(rep_constant){
      event_covs <- NULL
    } else {
      event_covs <- list(ec1 = t(utils::unstack(visit_full[c("ec1", "id_unit")], ec1 ~ id_unit))) 
      assertthat::assert_that(isTRUE(all.equal(rownames(obs), rownames(event_covs$ec1))))
    }

    ub2 <- unit_backbone[match(rownames(obs), paste0("X", unit_backbone$id_unit)), ]
    assertthat::assert_that(all.equal(rownames(obs), paste0("X", ub2$id_unit)))
    unit_covs = ub2[c("uc1", "species")]
    
    out <- list(obs = remove_rownames(obs), 
                unit_covs = remove_rownames(unit_covs),
                event_covs = lapply(event_covs, remove_rownames) # this yields an empty list if event_covs is NULL
                )
  } else if (augmented) {
    visit_full$id_unit <- interaction(visit_full$id_point, visit_full$species)
    unit_backbone$id_unit <- interaction(unit_backbone$id_point, unit_backbone$species)
    
    obs_temp <- list()
    counter <- 0
    for(i in 1:n_sp){
      ot <- t(utils::unstack(visit_full[visit_full$species == paste0("sp_", i), c("obs", "id_point")], obs ~ id_point))
      if (sum(ot, na.rm = TRUE) > 0) { # we only include species in the output if the species is observed
        counter <- counter + 1
        obs_temp[[counter]] <- ot
        assertthat::assert_that(all.equal(rownames(ot), rownames(obs_temp[[1]])))
      }

    }

    ec_prelim <- visit_full[visit_full$species == "sp_1", ]
    
    if(rep_constant){
      event_covs <- NULL
    } else {
      event_covs <- list(ec1 = t(utils::unstack(ec_prelim[c("ec1", "id_point")], ec1 ~ id_point))) 
      assertthat::assert_that(isTRUE(all.equal(rownames(obs_temp[[1]]), rownames(event_covs$ec1))))
    }
    
    # match just returns the positions of the first match
    ub2 <- unit_backbone[match(rownames(obs_temp[[1]]), paste0("X", unit_backbone$id_point)), ]
    assertthat::assert_that(all.equal(rownames(obs_temp[[1]]), paste0("X", ub2$id_point)))
    unit_covs = ub2[c("uc1")]
    
    out <- list(
      obs = remove_rownames(abind::abind(obs_temp, along = 3)),
      unit_covs = remove_rownames(unit_covs),
      event_covs = lapply(event_covs, remove_rownames)
    )
    
  } else if (n_season > 1) {
    obs_temp <- events_temp <- units_temp <- list()
    for(i in seq_len(n_season)){
      visit_full$id_unit <- interaction(visit_full$id_point, visit_full$species)
      unit_backbone$id_unit <- interaction(unit_backbone$id_point, unit_backbone$species)
      
      vfi <- visit_full[visit_full$id_season == i, ]
      
      obs_temp[[i]] <- t(utils::unstack(vfi[c("obs", "id_unit")], obs ~ id_unit))
      
      if(!rep_constant){
        events_temp[[i]] <- t(utils::unstack(vfi[c("ec1", "id_unit")], ec1 ~ id_unit))
        assertthat::assert_that(all.equal(rownames(obs_temp[[i]]), rownames(events_temp[[i]])))
      }
      
      assertthat::assert_that(all.equal(rownames(obs_temp[[1]]), rownames(obs_temp[[i]])))
      ub2 <- unit_backbone[match(rownames(obs_temp[[i]]), paste0("X", unit_backbone$id_unit)), ]
      assertthat::assert_that(all.equal(rownames(obs_temp[[i]]), paste0("X", ub2$id_unit)))
      
      units_temp[[i]] = ub2[c("uc1", "species")]
    }
    
    if(rep_constant){
      event_covs <- list()
    } else {
      event_covs <- list(ec1 = remove_rownames(abind::abind(events_temp, along = 3)))
    }
    
    
    out <- list(
      obs = remove_rownames(abind::abind(obs_temp, along = 3)),
      unit_covs = lapply(units_temp, remove_rownames),
      event_covs = event_covs
    )
  }

  out$params <- params
  out$covariates <- covariates
  out
}

