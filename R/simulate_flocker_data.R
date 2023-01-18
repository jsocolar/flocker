#' Simulate data for use with \code{make_flocker_data()} and downstream 
#' functions. 
#' 
#' Data will be simulated with one unit covariate that affects all 
#' relevant terms, one event covariate that affects detection (unless 
#' `rep_constant` is `TRUE`), and one grouping factor representing species
#' with correlated effects on all terms.
#' @param rep_constant logical: create data with unit covariates only (TRUE) 
#'   or data that includes event covariates (FALSE)
#' @param fp logical: if create data for fp model (TRUE) or standard model 
#'   (FALSE). If `TRUE`, the simulation will produce detections that correspond
#'   to true detections with probability 0.8.
#' @param n_pt number of points to simulate. The number of units for single-
#'   season models will be `n_pt*n_sp`. The number of units for multi-season
#'   models will be `n_pt*n_sp*n_season`.
#' @param n_sp number of levels to include in random effect. For compatibility
#'   with multispecies models, this levels of this random effect define separate
#'   closure units.
#' @param n_rep number of replicate visits to simulate per closure unit
#' @param ragged_rep logical: create data with variable (TRUE) or constant 
#'   (FALSE) numbers of visits per unit.  If TRUE, approximately half of units 
#'   will be missing approximately half of `n_rep` visits.
#' @param n_season A positive integer giving the number of seasons desired. 1
#'   yields data for a single-season model; all other positive integers yield
#'   multiseason models.
#' @param missing_seasons logical; relevant only if n_season is greater than 1. 
#'   create data with variable (TRUE) or constant (FALSE) numbers of seasons per 
#'   series (TRUE). If TRUE, approximately half of series will be missing their
#'   even-numbered seasons.
#' @param multiseason if n_season is NULL, must be NULL. Otherwise, one of
#'   "colex" or "autologistic"
#' @param multi_init if n_season is NULL, must be NULL. Otherwise, one of 
#'   "explicit" or "equilibrium"
#' @param augmented logical. If `TRUE` data will be formatted for an augmented
#'   model; `fp` must be `NULL`, and `n_season` must be `1`; all never-observed 
#'   species will be trimmed out of the data before augmenting; by default 
#'   data will be simulated under a large random effect variance for detection 
#'   to encourage the existence of never-observed species; and by default data 
#'   will not include any covariate effects on occupancy.
#' @param params a named list of parameter values to use in simulation, or NULL 
#'   to simulate values.
#' @param covariates a dataframe of covariate values to use in simulation, or
#'   NULL oto simulate values.
#' @param seed random seed. NULL uses (and updates) the existing RNG state.
#' @return A three element named list with the observation matrix/array ($obs), 
#'   the unit covariate dataframe(s) ($unit_covs), and the event covariate list
#'   ($rep_covs). If rep_constant is TRUE, then $rep_covs will be NULL.
#' @export

simulate_flocker_data <- function(
  n_rep = 4,
  n_pt = 50,
  n_sp = 30,
  n_season = 1,
  multiseason = NULL,
  multi_init = NULL,
  augmented = FALSE,
  rep_constant = FALSE,
  fp = NULL,
  params = NULL,
  covariates = NULL,
  seed = 123,
  ragged_rep = FALSE,
  missing_seasons = FALSE) {
  if (is.null(seed)){
    out <- sfd(n_rep, n_pt, n_sp, n_season, multiseason, multi_init, augmented,
               rep_constant, fp, params, covariates, ragged_rep, missing_seasons)
  } else {
      out <- withr::with_seed(
        seed,
        sfd(n_rep, n_pt, n_sp, n_season, multiseason, multi_init, augmented,
            rep_constant, fp, params, covariates, ragged_rep, missing_seasons)
      )
  }
  out
}
  
#' util for creating example data
#' @inheritParams simulate_flocker_data
#' @return A three element named list with the observation matrix/array ($obs), 
#'   the unit covariate dataframe(s) ($unit_covs), and the event covariate list
#'   ($rep_covs). If rep_constant is TRUE, then $rep_covs will be NULL.
sfd <- function(
    n_rep,
    n_pt,
    n_sp,
    n_season,
    multiseason,
    multi_init,
    augmented,
    rep_constant,
    fp,
    params,
    covariates,
    ragged_rep,
    missing_seasons
) {
  assertthat::assert_that(
    n_rep == floor(n_rep) & n_rep > 1,
    msg = "n_rep must be an integer greater than 1"
    )
  
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
        Sigma[which(coef_names == "det_intercept"), ] <- 3 * Sigma[which(coef_names == "det_intercept"), ]
        Sigma[, which(coef_names == "det_intercept")] <- 3 * Sigma[, which(coef_names == "det_intercept")]
        Sigma[which(coef_names == "occ_intercept"), ] <- 3 * Sigma[which(coef_names == "occ_intercept"), ]
        Sigma[, which(coef_names == "occ_intercept")] <- 3 * Sigma[, which(coef_names == "occ_intercept")]
        Sigma[which(coef_names == "det_slope_unit"), ] <- .1 * Sigma[which(coef_names == "det_slope_unit"), ]
        Sigma[, which(coef_names == "det_slope_unit")] <- .1 * Sigma[, which(coef_names == "det_slope_unit")]
        if(!rep_constant) {
          Sigma[which(coef_names == "det_slope_visit"), ] <- .1 * Sigma[which(coef_names == "det_slope_visit"), ]
          Sigma[, which(coef_names == "det_slope_visit")] <- .1 * Sigma[, which(coef_names == "det_slope_visit")]
        }
      }
      params$Sigma <- Sigma
    }
    params$coefs <- as.data.frame(
      MASS::mvrnorm(n_sp, params$coef_means, params$Sigma)
    )
    names(params$coefs) <- coef_names
  }
  assertthat::assert_that(
    identical(names(params$coefs), coef_names)
  )
  coefs <- params$coefs
  if (isTRUE(multiseason == "autologistic")) {
    if (!("theta" %in% np)) {
      params$theta <- 1
    }
  }
  
  # Create covariates (currently assume that unit covariates are the same
  # across all timesteps in multiseason models)
  unit_backbone <- expand.grid(
    species = factor(paste0("sp_", 1:n_sp), levels = paste0("sp_", 1:n_sp)),
    id_point = 1:n_pt
    )
  if (is.null(covariates)) {
    unit_backbone$uc1 <- rnorm(n_pt)[unit_backbone$id_point]
  } else {
    unit_backbone$uc1 <- covariates$uc1
  }
  
  if (n_season > 1) {
    unit_backbone$col <- coefs$col_intercept[as.integer(unit_backbone$species)] +
      coefs$col_slope_unit[as.integer(unit_backbone$species)] * unit_backbone$uc1
    if (multiseason == "colex") {
      unit_backbone$ex <- coefs$ex_intercept[as.integer(unit_backbone$species)] +
        coefs$ex_slope_unit[as.integer(unit_backbone$species)] * unit_backbone$uc1
    } else if (multiseason == "autologistic") {
      unit_backbone$ex <- boot::logit(
        1 - boot::inv.logit(params$theta + unit_backbone$col)
      )
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
      true_Z = rbinom(n_sp * n_pt, 1, boot::inv.logit(unit_backbone$occ))
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
              true_Z = rbinom(n_pt * n_sp, 1, psi_i)
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
  
  visit_full <- merge(
    visit_backbone, 
    true_Z, 
    by = c("species", "id_point", "id_season"),
    all = TRUE
  )
  
  visit_full$logit_det <- 
    coefs$det_intercept[as.integer(visit_full$species)] +
    coefs$det_slope_unit[as.integer(visit_full$species)] * visit_full$uc1
  
  if (!rep_constant) {
    if (is.null(covariates)) {
      visit_full$vc1 <- rnorm(n_pt*n_season*n_rep)[ # we don't allow vc1 to vary by species;
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
    } else {
      visit_full$vc1 <- covariates$vc1
    }
    visit_full$logit_det <- visit_full$logit_det + 
      coefs$det_slope_visit[as.integer(visit_full$species)] * visit_full$vc1
  }
  
  visit_full$obs <- visit_full$true_Z *
    rbinom(nrow(visit_full), 1, boot::inv.logit(visit_full$logit_det))
  
  if (!is.null(fp)) {
    # true detections are categorized as such with probability fp
    n_succ <- sum(visit_full$obs)
    visit_full$obs[visit_full$obs == 1] <- fp
    # other detections are categorized as {true detections with prob. fp} with
    # probability so as to imply that on average a fraction fp of detections 
    # are true
    n_fail <- stats::rnbinom(1, n_succ, fp) # this is approximate, but a pretty good one
    n_zeros <- sum(!visit_full$obs) # 
    if(n_fail > n_zeros) {
      stop(
      "there are not enough true nondetections to accommodate the desired fp probability"
      )
    }
    zeros <- rep(0, n_zeros)
    zeros[sample(1:n_zeros, n_fail)] <- fp
    visit_full$obs[visit_full$obs == 0] <- zeros
  }
  
  if (ragged_rep) {
    visit_full$season_pt <- interaction(visit_full$id_point, visit_full$id_season)
    season_pt_missing <- sample(unique(visit_full$season_pt), floor(length(unique(visit_full$season_pt)) / 2))
    rows_missing <- which((visit_full$season_pt %in% season_pt_missing) & (visit_full$id_rep > n_rep / 2))
    visit_full$obs[rows_missing] <- visit_full$vc1[rows_missing] <- NA
  }
  
  if (missing_seasons) {
    pt_missing <- sample(unique(visit_full$id_point), floor(length(unique(visit_full$id_point)) / 2))
    rows_missing <- which((visit_full$season_pt %in% pt_missing) & (visit_full$id_season %% 2 == 0))
    visit_full$obs[rows_missing] <- visit_full$vc1[rows_missing] <- NA
  }
  
  # format data for return
  if (! augmented & n_season == 1) {
    visit_full$id_unit <- interaction(visit_full$id_point, visit_full$species)
    unit_backbone$id_unit <- interaction(unit_backbone$id_point, unit_backbone$species)
    obs <- t(unstack(visit_full[c("obs", "id_unit")], obs ~ id_unit))
    event_covs = list(
      vc1 = t(unstack(visit_full[c("vc1", "id_unit")], vc1 ~ id_unit))
    )
    assertthat::assert_that(all.equal(rownames(obs), rownames(event_covs$vc1)))
    ub2 <- unit_backbone[match(rownames(obs), paste0("X", unit_backbone$id_unit)), ]
    assertthat::assert_that(all.equal(rownames(obs), paste0("X", ub2$id_unit)))
    
    unit_covs = ub2[c("uc1", "species")]
    
    out <- list(obs = remove_rownames(obs), 
                unit_covs = remove_rownames(unit_covs),
                event_covs = lapply(event_covs, remove_rownames)
                )
  } else if (augmented) {
    visit_full$id_unit <- interaction(visit_full$id_point, visit_full$species)
    unit_backbone$id_unit <- interaction(unit_backbone$id_point, unit_backbone$species)
    
    obs_temp <- list()
    counter <- 0
    for(i in 1:n_sp){
      ot <- t(unstack(visit_full[visit_full$species == paste0("sp_", i), c("obs", "id_point")], obs ~ id_point))
      if (sum(ot) > 0) {
        counter <- counter + 1
        obs_temp[[counter]] <- ot
        assertthat::assert_that(all.equal(rownames(ot), rownames(obs_temp[[1]])))
      }

    }

    ec_prelim <- visit_full[visit_full$species == "sp_1", ]
    event_covs = list(
      vc1 = t(unstack(ec_prelim[c("vc1", "id_point")], vc1 ~ id_point))
    )
    assertthat::assert_that(all.equal(rownames(obs_temp[[1]]), rownames(event_covs$vc1)))
    
    
    ub2 <- unit_backbone[match(rownames(obs_temp[[1]]), paste0("X", unit_backbone$id_point)), ]
    assertthat::assert_that(all.equal(rownames(obs_temp[[1]]), paste0("X", ub2$id_point)))
    
    unit_covs = ub2[c("uc1", "species")]
    
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
      
      obs_temp[[i]] <- t(unstack(vfi[c("obs", "id_unit")], obs ~ id_unit))
      
      events_temp[[i]] <- t(unstack(vfi[c("vc1", "id_unit")], vc1 ~ id_unit))
      
      assertthat::assert_that(all.equal(rownames(obs_temp[[1]]), rownames(obs_temp[[i]])))
      assertthat::assert_that(all.equal(rownames(obs_temp[[i]]), rownames(events_temp[[i]])))
      ub2 <- unit_backbone[match(rownames(obs_temp[[i]]), paste0("X", unit_backbone$id_unit)), ]
      assertthat::assert_that(all.equal(rownames(obs_temp[[i]]), paste0("X", ub2$id_unit)))
      
      units_temp[[i]] = ub2[c("uc1", "species")]
    }
    out <- list(
      obs = remove_rownames(abind::abind(obs_temp, along = 3)),
      unit_covs = lapply(units_temp, remove_rownames),
      event_covs = list(vc1 = remove_rownames(abind::abind(events_temp, along = 3)))
      )
  }

  out
}

