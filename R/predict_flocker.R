#' Get posterior predictions from a flocker model
#' @param flocker_fit A `flocker_fit` object
#' @param draw_ids Vector of indices of the posterior draws to be 
#'     used. If `NULL` (the default) all draws are used in their native order.
#' @param history_condition Logical indicator of whether to directly condition the 
#'     posterior for the occupancy state on the observed detection histories.
#'     For example, at sites with at least one detection, the true occupancy 
#'     state conditioned on the history is one with absolute certainty. Without 
#'     directly conditioning on the history, the occupancy state is controlled 
#'     exclusively by the posterior distribution for the occupancy probability 
#'     psi.
#' @param new_data Optional new data at which to predict. If `NULL`, predictions
#'     are given at the data points used for model fitting. Otherwise, must
#'     be a flocker_data object produced by `make_flocker_data`.
#' @param mixed When `new_data` is not provided, should random effect levels be
#'     drawn from their posteriors (`FALSE`, the default) or re-sampled from 
#'     their fitted hyperparameters (`TRUE`). The latter can be useful for mixed
#'     predictive checking. Note that setting to TRUE is not available for
#'     grouping terms involved in phylogenetic random effects or other random 
#'     effects with specified covariance structures.
#' @param allow_new_levels Should new_data be allowed to contain new levels for
#'     random effects?
#' @param sample_new_levels If new_data is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to brms::prepare_predictions, which see. 
#' @return An array of posterior predictions in the same shape as the
#'     observations passed to `make_flocker_data()` with posterior iterations
#'     stacked along the final dimension.
#' @export
#' @examples 
#' \donttest{
#' unconditioned_preds <- predict_flocker(example_flocker_model_single)
#' conditioned_preds <- predict_flocker(
#'  example_flocker_model_single, 
#'  history_condition = TRUE
#' )
#' }
predict_flocker <- function(flocker_fit, draw_ids = NULL,  
                            history_condition = FALSE, 
                            new_data = NULL,
                            mixed = FALSE, 
                            allow_new_levels = FALSE,
                            sample_new_levels = "uncertainty") {
  assertthat::assert_that(
    is_one_logical(mixed),
    msg = "`mixed` must be a single logical"
  )
  
  assertthat::assert_that(is.null(new_data) | is_flocker_data(new_data))
  
  total_iter <- brms::ndraws(flocker_fit)
  
  if (is.null(draw_ids)) {
    n_iter <- total_iter
  } else {
    n_iter <- length(draw_ids)
  }
  
  # rename all random effect levels so they show up as new levels
  if (mixed) {
    if (is.null(new_data)) {
      new_data <- flocker_fit$data
    } else {
      assertthat::assert_that(
        identical(new_data$type, attributes(flocker_fit)$data_type),
        msg = "the new_data data type does not match the flocker_fit data type"
      )
      new_data <- new_data$data
    }
    
    random_effects <- unique(flocker_fit$ranef$group)
    if (length(random_effects) > 0) {
      potential_conflicts <- vector()
      for(i in 1:length(random_effects)) {
        potential_conflicts <- c(potential_conflicts, unique(paste0(random_effects[i], flocker_fit$data[, random_effects[i]])))
      }
      potential_conflicts <- c(potential_conflicts, paste0("occ_", potential_conflicts))
      fixed_effects <- rownames(brms::fixef(flocker_fit))
      if (any(potential_conflicts %in% fixed_effects)) {
        stop(paste0("Mixed predictions not implemented when some groups ",
                    "appear as both fixed and random groupings in different parts ",
                    "of the model. Posterior predictions with `mixed = FALSE` ",
                    "are enabled for this model."))
      }
      for (i in 1:length(random_effects)) {
        new_data[, random_effects[i]] <- paste0(new_data[, random_effects[i]], 
                                              "_resampled")
      }
      new_data <- list(data = new_data, type = attributes(flocker_fit)$data_type)
      class(new_data) <- "flocker_data"
    }
    assertthat::assert_that(
      identical(allow_new_levels, TRUE),
      msg = "set `allow_new_levels` to TRUE for mixed predictive checking"
    )
    assertthat::assert_that(
      identical(sample_new_levels, "gaussian"),
      msg = "set `sample_new_levels` to 'gaussian' for mixed predictive checking"
      )
  }
  
  
  Z_samp <- get_Z(flocker_fit, draw_ids = draw_ids, history_condition = history_condition, 
                  sample = TRUE, new_data = new_data, 
                  allow_new_levels = allow_new_levels, sample_new_levels = sample_new_levels)
  
  lps <- fitted_flocker(
    flocker_fit, 
    components = "det",
    draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
    sample_new_levels = sample_new_levels, response = FALSE, unit_level = FALSE
  )
  theta_all <- boot::inv.logit(lps$linpred_det)
  ndim <- length(dim(theta_all))
  
  assertthat::assert_that(
    ndim > 2,
    msg = "predict_flocker error 1. This shouldn't happen; please report a bug"
  )
  
  Z_samp_array <- abind::abind(rep(list(Z_samp), dim(theta_all)[2]), along = ndim) |>
    aperm(perm = c(1, ndim, (2:(ndim - 1))))
  
  assertthat::assert_that(
    identical(dim(theta_all), dim(Z_samp_array)),
    msg = "predict_flocker error 2. This shouldn't happen; please report a bug"
  )
  
  predictions <- new_array(theta_all, rbinom2(length(theta_all), 1, theta_all * Z_samp_array))
  
  predictions
}
