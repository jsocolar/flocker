#' Get posterior predictions from a flocker model
#' @param flocker_fit A `flocker_fit` object
#' @param n_iter The number of posterior iterations desired. If `NULL`, use
#'     all available posterior iterations.
#' @param new_data Optional new data at which to predict. If `NULL`, predictions
#'     are given at the data points used for model fitting ("retrodictions")
#' @param history_condition Logical indicator of whether to directly condition the 
#'     posterior for the occupancy state on the observed detection histories.
#'     For example, at sites with at least one detection, the true occupancy 
#'     state conditioned on the history is one with absolute certainty. Without 
#'     directly conditioning on the history, the occupancy state is controlled 
#'     by the posterior distribution for the occupancy probability psi. Of 
#'     course even without conditioning directly on the detection history, we 
#'     still condition indirectly on the observed history via the fitted value 
#'     of psi, which itself depends on all of the observed detection histories.
#' @param mixed When `new_data` is not provided, should random effect levels be
#'     drawn from their posteriors (`FALSE`, the default) or re-sampled from 
#'     their fitted hyperparameters (`TRUE`). The latter can be useful for mixed
#'     predictive checking.
#' @param sample_new_levels If new_data is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to brms::prepare_predictions, which see. 
#' @return An array of posterior predictions. If the model is rep-varying, 
#'     then a 3-dimensional array where the first dimension is the closure-unit,
#'     the second dimension is the rep, the third dimension is the iteration, 
#'     and the value is 1, 0, or NA indicating detection, non-detection, or 
#'     non-existence of the sampling event.
#'     If the model is rep-constant, then a matrix where rows are iterations,
#'     columns are closure-units, and values are the number of successes.
#' @export

predict_flocker <- function(flocker_fit, n_iter=NULL, new_data = NULL, 
                            history_condition = FALSE, mixed = FALSE, 
                            sample_new_levels = "uncertainty") {
  if (!is_flocker_fit(flocker_fit)) {
    stop("`flocker_fit` must be an object of class `flocker_fit`")
  }
  if (!is.logical(history_condition)) {
    stop("`history_condition` must be logical")
  }
  if (history_condition & (!is.null(new_data))) {
    stop("`history_condition` must be `FALSE` if new_data is supplied")
  }
  if (!is.logical(mixed)) {
    stop("`mixed` must be logical")
  }
  if (mixed & (!is.null(new_data))) {
    stop("`mixed` must be `FALSE` if new_data is supplied")
  }
  new_data2 <- new_data
  
  if (is.null(new_data)) {
    new_data <- flocker_fit$data
  } else {
    new_data <- new_data$data
  }
  
  total_iter <- brms::niterations(flocker_fit)*brms::nchains(flocker_fit)
  
  if (is.null(n_iter)) {
    n_iter <- total_iter
  }
  
  n_iter <- as.integer(n_iter)
  
  if ((length(n_iter) != 1) | is.na(n_iter)) {
    stop("iter, if supplied, must be a single integer")
  }
  
  if (n_iter > total_iter) {
    stop("requested more iterations than contained in flocker_fit")
  }
  
  if (mixed) {
    random_effects <- flocker_fit$ranef$group
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
    }
    sample_new_levels = "gaussian"
    message("`sample_new_levels` set to 'gaussian' for mixed predictive checking")
  }
  
  Z <- get_Z(flocker_fit, n_iter = n_iter, history_condition = history_condition, 
             new_data = new_data2, sample_new_levels = sample_new_levels)
  Z_samp <- apply(Z, 2, function(x){rbinom(length(x), 1, x)})
  
  iter <- (1:n_iter)*floor(total_iter/n_iter)
  
  # The below line is inefficient because it duplicates work already performed
  # by `get_Z`, above. A potential target for refactoring.
  pd <- boot::inv.logit(
    brms::posterior_linpred(flocker_fit, dpar = "mu", draw_ids = iter, 
                            newdata = new_data, allow_new_levels = TRUE, 
                            sample_new_levels = sample_new_levels)
  )
  
  lt <- type_flocker_fit(flocker_fit)
  
  if (lt == "single") {
    n_unit <- new_data$n_unit[1]
    max_rep <- max(new_data$n_rep)
    positions <- get_positions_V(flocker_fit)
    out <- array(dim = c(n_unit, max_rep, n_iter), 
                 dimnames = list(unit = NULL, rep = NULL, iter = NULL))
    for (k in 1:nrow(flocker_fit$data)) {
      i <- positions[k, 1]
      j <- positions[k, 2]
      out[i, j, ] <- rbinom(n_iter, 1, pd[, k]*Z_samp[, i])
    }
  } else if (lt == "single_C") {
    n_unit <- nrow(new_data)
    out <- matrix(nrow = n_iter, ncol = n_unit)
    for (i in 1:n_unit) {
      out[, i] <- rbinom(n_iter, flocker_fit$data$n_trial[i], pd[, i]*Z_samp[, i])
    }
  }
  return(out)
}
