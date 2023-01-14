#' Create example data for use with \code{make_flocker_data()} and downstream 
#' functions. Data will be simulated with one unit covariate that effects all 
#' relevant terms, one event covariate that affects detection (unless 
#' `rep_constant` is `TRUE`), and one grouping factor with correlated effects on 
#' all terms.
#' @param rep_constant logical: create data with unit covariates only (TRUE) 
#'   or data that includes event covariates (FALSE)
#' @param fp logical: if create data for fp model (TRUE) or standard model 
#'   (FALSE). If `TRUE`, all simulated detections will be relabeled as having 
#'   a priori probability of 0.8 of corresponding to true detections.
#' @param n_pt number of points to simulate.
#' @param n_sp number of levels to include in random effect, which will be
#'   treated as "species" if an augmented model is requested via `omega`.
#'   For an augmented model, this is the expected number of species in the
#'   community, but is not guaranteed to be the exact number.
#' @param n_rep number of replicate visits to simulate per closure unit
#' @param ragged_rep logical: create data with variable (TRUE) or constant 
#'   (FALSE) numbers of visits per unit.  If TRUE, approximately half of units 
#'   will be missing approximately half of `n_rep` visits.
#' @param n_season if NULL, data for a single season model are returned.
#'   Otherwise, an integer giving the number of seasons desired in for a 
#'   multiseason model.
#' @param missing_seasons logical; relevant only if n_season is specified: 
#'   create data with variable (TRUE) or constant (FALSE) numbers of seasons per 
#'   series (TRUE). If TRUE, approximately half of series will be missing their
#'   even-numbered seasons.
#' @param multi_type if n_season is NULL, must be NULL. Otherwise, one of
#'   "colex", "equilibrium"
#' @param omega NULL for models other than the data-augmented model. For the
#'   data augmented model, a number between 0 and 1 giving the probability 
#'   that a species from the augmented set is included in the community. 
#'   Must be NULL if either `fp` or `n_season` is not NULL. If TRUE, by default
#'   data will be simulated under a large random effect variance for
#'   detection, which should encourage the existence of some never-observed
#'   species; and by default data will not include any unit covariate effects
#'   on occupancy.
#' @param params a list of parameter values to use in simulation, or NULL to
#'   use default values. If not NULL, it is the user's responsibility to ensure
#'   that values are passed for all relevant parameters.
#' @param seed random seed. NULL uses (and updates) the existing RNG state.
#' @return A three element named list with the observation matrix/array ($obs), 
#'   the unit covariate dataframe(s) ($unit_covs), and the event covariate list
#'   ($rep_covs). If rep_constant is TRUE, then $rep_covs will be NULL.
#' @export

example_flocker_data <- function(
  rep_constant = FALSE,
  fp = NULL,
  n_pt = 30,
  n_sp = 30,
  n_rep = 4,
  ragged_rep = FALSE,
  n_season = NULL,
  missing_seasons = FALSE,
  augmented = FALSE,
  params = NULL,
  seed = NULL) {
  if (is.null(seed)){
    out <- efd(rep_constant, fp, n_pt, n_sp, n_rep, ragged_rep, n_season,
      missing_seasons, n_aug, omega, verbose)
  } else {
      out < withr::with_seed(
        seed,
        efd(rep_constant, fp, n_pt, n_sp, n_rep, ragged_rep, n_season,
          missing_seasons, n_aug, omega, verbose)
      )
  }
  out
}
  
#' util for creating example data
#' @inheritParams example_flocker_data
#' @return A three element named list with the observation matrix/array ($obs), 
#'   the unit covariate dataframe(s) ($unit_covs), and the event covariate list
#'   ($rep_covs). If rep_constant is TRUE, then $rep_covs will be NULL.
efd <- function(
    rep_constant,
    fp,
    n_pt,
    n_sp,
    n_rep,
    ragged_rep,
    n_season,
    missing_seasons,
    augmented,
    params
) {
  n_coef <- 3 + !rep_constant +
  if(is.null(params)) {
  }
  Sigma <- matrix(.5, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  if (rep_constant) {
    
    
  }
    
    
  MASS::mvrnorm(n_sp, rep(0, 4), matrix(c(1,.5,.5,.5,
                                       .5,1,.5,.5,
                                       .5,.5,1,.5,
                                       .5,.5,.5,1), nrow = 4))
}





  n_unit <- n_pt*n_sp
  backbone <- expand.grid(species = factor(paste0("sp_", 1:n_sp)), 
                          id_rep = 1:n_rep,
                          id_point = 1:n_pt)
  
  unit_covs <- within(unique(backbone[c("species", "id_point")]), {
    id_unit = 1:n_unit
    uc1 = stats::rnorm(n_pt)[id_point]
    uc2 = stats::rnorm(n_pt)[id_point]
    grp = sample(c(1:20), n_unit, replace = T)
    logit_psi = 0 + 1*uc1 - 1*uc2 + stats::rnorm(20)[grp]
    true_Z = stats::rbinom(n_unit, 1, boot::inv.logit(logit_psi))
  })
  
  df_full <- within(merge(backbone, unit_covs), {
    id_point_rep <- interaction(id_point, id_rep)
    ec1 = stats::rnorm(n_pt * n_rep)[id_point_rep]
    ec2 = stats::rnorm(n_pt * n_rep)[id_point_rep]
    logit_theta = -.5 + .5*uc1 + stats::rnorm(20)[grp]
    if (!rep_constant) {
      logit_theta = logit_theta + .5*ec1 - 1*ec2
    }
    obs <- stats::rbinom(n_unit * n_rep, true_Z, boot::inv.logit(logit_theta))
  })
  
  # format data for return
  out <- list(obs = t(unstack(df_full[c("obs", "id_unit")], obs ~ id_unit)), 
              unit_covs = unit_covs[c("uc1", "uc2", "grp", "species")], 
              event_covs = list(
                ec1 = t(unstack(df_full[c("ec1", "id_unit")], ec1 ~ id_unit)), 
                ec2 = t(unstack(df_full[c("ec2", "id_unit")], ec2 ~ id_unit))
              ))
  rownames(out$obs) <- NULL
  rownames(out$event_covs$ec1) <- NULL
  rownames(out$event_covs$ec2) <- NULL
  return(out)

