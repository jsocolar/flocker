##### single #####

#' Create Stan code for likelihood function occupancy_single_lpmf.
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @return Character string of Stan code corresponding to occupancy_single_lpmf
#' @noRd
make_occupancy_single_lpmf <- function (max_rep) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_single_lpmf(
    array[] int y, // detection data
    vector mu, // lin pred for detection
    vector occ, // lin pred for occupancy. Elements after vint1[1] irrelevant.
    array[] int vint1, // # units (n_unit). Elements after 1 irrelevant.
    array[] int vint2, // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
    array[] int vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
  
  // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):"
  
  sf_text2 <- paste0("    array[] int vint", 3 + (1:max_rep), collapse = ",\n")
  
  sf_text3 <- paste0(") {
  // Create array of the rep indices that correspond to each unit.
    array[vint1[1], ", max_rep, "] int index_array;")
  
  sf_text4.1 <- "      index_array[,"
  sf_text4.2 <- 1:max_rep
  sf_text4.3 <- "] = vint"
  sf_text4.4 <- 3 + (1:max_rep)
  sf_text4.5 <- "[1:vint1[1]];\n"
  sf_text4 <- paste0(sf_text4.1, sf_text4.2, sf_text4.3, sf_text4.4, sf_text4.5, collapse = "")
  
  sf_text5 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    for (i in 1:vint1[1]) {
      array[vint2[i]] int indices = index_array[i, 1:vint2[i]];
      if (vint3[i] == 1) {
        lp += bernoulli_logit_lpmf(1 | occ[i]);
        lp += bernoulli_logit_lpmf(y[indices] | mu[indices]);
      }
      if (vint3[i] == 0) {
        lp += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                              sum(log1m_inv_logit(mu[indices])), bernoulli_logit_lpmf(0 | occ[i]));
      }
    }
    return(lp);
  }
"

  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sep = "\n")
  out
}

##### single_C #####

#' Create Stan code for likelihood function occupancy_single_C_lpmf.
#' The purpose of defining this custom family, rather than using brms's zero-inflated 
#' binomial, is to ensure that the occupancy parameters are interpretable with
#' values of 1 in the marginalized state reflecting occupancy rather than non-
#' occupancy.
#' @return Character string of Stan code corresponding to occupancy_single_C_lpmf
#' @noRd
make_occupancy_single_C_lpmf <- function () {
  "real occupancy_single_C_lpmf(int y, real mu, real occ, int trials) {
  if (y == 0) { 
    return log_sum_exp(bernoulli_logit_lpmf(0 | occ), 
                       bernoulli_logit_lpmf(1 | occ) + 
                         binomial_logit_lpmf(0 | trials, mu)); 
  } else { 
    return bernoulli_logit_lpmf(1 | occ) +  
      binomial_logit_lpmf(y | trials, mu); 
  } 
}"
}

##### augmented #####

#' Create Stan code for likelihood function occupancy_augmented_lpmf for 
#' rep-varying model. 
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @return Character string of Stan code corresponding to occupancy_augmented_lpmf
#' @noRd
make_occupancy_augmented_lpmf <- function (max_rep) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_augmented_lpmf(
    array[] int y, // detection data
    vector mu, // lin pred for detection
    vector occ, // lin pred for occupancy. Elements after vint1[1] irrelevant.
    vector Omega, // lin pred for availability.  Elements after 1 irrelevant.
    array[] int vint1, // # units (n_unit). Elements after 1 irrelevant.
    array[] int vint2, // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
    array[] int vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
    
    array[] int vint4, // # species (observed + augmented). Elements after 1 irrelevant.
    array[] int vint5, // Indicator for species was observed.  Elements after vint4[1] irrelevant
    
    array[] int vint6, // species
  
  // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):"
  
  sf_text2 <- paste0("    array[] int vint", 6 + (1:max_rep), collapse = ",\n")
  
  sf_text3 <- paste0(") {
  // Create array of the rep indices that correspond to each unit.
    array[vint1[1], ", max_rep, "] int index_array;")
  
  sf_text4.1 <- "      index_array[,"
  sf_text4.2 <- 1:max_rep
  sf_text4.3 <- "] = vint"
  sf_text4.4 <- 6 + (1:max_rep)
  sf_text4.5 <- "[1:vint1[1]];\n"
  sf_text4 <- paste0(sf_text4.1, sf_text4.2, sf_text4.3, sf_text4.4, sf_text4.5, collapse = "")
  
  sf_text5 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    
    for (sp in 1:vint4[1]) {
      real lp_s = 0;
      if (vint5[sp] == 1) {
        for (i in 1:vint1[1]) {
          if (vint6[i] == sp) {
            array[vint2[i]] int indices = index_array[i, 1:vint2[i]];
            if (vint3[i] == 1) {
              lp_s += bernoulli_logit_lpmf(1 | occ[i]);
              lp_s += bernoulli_logit_lpmf(y[indices] | mu[indices]);
            }
            if (vint3[i] == 0) {
              lp_s += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                                    sum(log1m_inv_logit(mu[indices])), bernoulli_logit_lpmf(0 | occ[i]));
            }
          }
        }
        lp += log_inv_logit(Omega[1]) + lp_s;
      } else {
        for (i in 1:vint1[1]) {
          if (vint6[i] == sp) {
            array[vint2[i]] int indices = index_array[i, 1:vint2[i]];
            lp_s += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                                  sum(log1m_inv_logit(mu[indices])), bernoulli_logit_lpmf(0 | occ[i]));
          }
        }
        lp += log_sum_exp(log1m_inv_logit(Omega[1]), log_inv_logit(Omega[1]) + lp_s);  
      }
    }
    return(lp);
  }
"
  
  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sep = "\n")
  return(out)
}


##### emission probabilities #####
# These are used primarily in multiseason models

# We don't write out the emission likelihood for zeros 
# because we use only the trivial likelihood associated with a 
# zero in the data and we drop mismatched terms explicitly in the way we
# compute the forward algorithm

#' Create Stan code for the emission log-likelihood given that the true state 
#'   equals one in a rep-varying model.
#' @return Character string of Stan code corresponding to emission_1
#' @noRd
make_emission_1 <- function() {
  paste("  // emission likelihood given that state equals one",
        "  real emission_1(array[] int y, row_vector det) {",
        "    return(bernoulli_logit_lpmf(y | det));",
        "  }",
        sep = "\n")
}

#' Create Stan code for the emission log-likelihood given that the true state 
#'   equals one in a rep-constant model.
#' @return Character string of Stan code corresponding to emission_1
#' @noRd
make_emission_1_C <- function() {
  paste("  // emission likelihood given that state equals one",
        "  real emission_1(int y, int k, real det) {",
        "    return(binomial_logit_lpmf(y | k, det));",
        "  }",
        sep = "\n")
}

##### multi-season: transition probabilities #####

#' Create Stan code for transition log-likelihoods in colex model.
#' These likelihoods ignore disallowed possibilities (detections on unoccupied
#' units); such possibilities are excluded via the control flow in the 
#' implementation of the forward algorithm.
#' @return Character string of Stan code corresponding to likelihoods for the 
#'  possible transitions in a colex model
#' @noRd
make_colex_likelihoods <- function() {
  paste("  // functions for the yearwise likelihood components corresponding to each",
        "  // possible pair of true states in that year and the previous year.",
        "  // These likelihoods drop the disallowed possibilities (detections",
        "  // on non-occupied units), which are excluded via the control flow",
        "  // in forward_colex.",
        "  real zero_zero(real colo) {",
        "    return(bernoulli_logit_lpmf(0 | colo));",
        "  }",
        "  real one_zero(real ex) {",
        "    return(bernoulli_logit_lpmf(1 | ex));",
        "  }",
        "  real zero_one(real colo, real e1) {  // e1 is the emission likelihood conditional on state == 1",
        "    return(bernoulli_logit_lpmf(1 | colo) + e1);",
        "  }",
        "  real one_one(real ex, real e1) {",
        "    return(bernoulli_logit_lpmf(0 | ex) + e1);",
        "  }",
        sep = "\n"
  )
}

##### multi-season: forward algorithm #####

#' Create Stan code for forward_colex function, an implementation of the 
#'  forward algorithm for use in colex models.
#' @return Character string of Stan code corresponding to forward_colex
#' @noRd
make_forward_colex <- function() {
  paste(
    "  // Forward algorithm implementation",
    "  real forward_colex(",
    "    int n_year, // number of years",
    "    array[] int Q, // At least one detection in year?",
    "    array[] int n_obs, // number of visits in year",
    "    array[,] int y, // detection nondetection data: rows are years, columns visits",
    "    real occ_initial, // initial occupancy logit-probability",
    "    vector colo, // colonization logit-probabilities (with initial dummy)",
    "    vector ex, // extinction logit-probabilities (with initial dummy)",
    "    array[] row_vector det // detection logit-probabilities",
    "  ) {",
    "    real alpha_0;",
    "    real alpha_1;",
    "    real alpha_0_prev;",
    "    real alpha_1_prev;",
    "    real e1; // emission probability conditional on latent state being 1.",
    "  ",
    "    // Year 1",
    "    if (n_obs[1] == 0) {",
    "      alpha_0 = bernoulli_logit_lpmf(0 | occ_initial);",
    "      alpha_1 = bernoulli_logit_lpmf(1 | occ_initial);",
    "    } else {",
    "      e1 = emission_1(y[1, 1:n_obs[1]], det[1, 1:n_obs[1]]);",
    "      if (Q[1] == 0) {",
    "        alpha_0 = bernoulli_logit_lpmf(0 | occ_initial);",
    "      }",
    "      alpha_1 = bernoulli_logit_lpmf(1 | occ_initial) + e1;",
    "    }",
    "    ",
    "    // Recursion for subsequent years",
    "    if (n_year > 1) {",
    "      // alpha_0  is not computed if Q[i] == 1",
    "      // alpha_1 is always computed.",
    "      for (i in 2:n_year) {",
    "        // Store alpha_0 and alpha_1 for the next round",
    "        alpha_0_prev = alpha_0;",
    "        alpha_1_prev = alpha_1;",
    "        ",
    "        if (n_obs[i] == 0) {  // year with no observations",
    "          if (Q[i - 1] == 0) {",
    "            alpha_0 = log_sum_exp(alpha_0_prev + bernoulli_logit_lpmf(0 | colo[i]),",
    "                                  alpha_1_prev + bernoulli_logit_lpmf(1 | ex[i]));",
    "            alpha_1 = log_sum_exp(alpha_0_prev + bernoulli_logit_lpmf(1 | colo[i]),",
    "                                  alpha_1_prev + bernoulli_logit_lpmf(0 | ex[i]));",
    "          } else if (Q[i - 1] == 1) {",
    "            alpha_0 = alpha_1_prev + bernoulli_logit_lpmf(1 | ex[i]);",
    "            alpha_1 = alpha_1_prev + bernoulli_logit_lpmf(0 | ex[i]);",
    "          }",
    "           ",  
    "        } else if (n_obs[i] > 0) { // year with observations",
    "          e1 = emission_1(y[i, 1:n_obs[i]], det[i, 1:n_obs[i]]);",
    "          if ((Q[i - 1] == 0) && (Q[i] == 0)) { // now years with at least one observation",
    "            alpha_0 = log_sum_exp(alpha_0_prev + zero_zero(colo[i]),",
    "                                  alpha_1_prev + one_zero(ex[i]));",
    "                                  alpha_1 = log_sum_exp(alpha_0_prev + zero_one(colo[i], e1),",
    "                                                        alpha_1_prev + one_one(ex[i], e1));",
    "          } else if ((Q[i - 1] == 0) && (Q[i] == 1)) {",
    "            alpha_1 = log_sum_exp(alpha_0_prev + zero_one(colo[i], e1),",
    "                                  alpha_1_prev + one_one(ex[i], e1));",
    "          } else if ((Q[i - 1] == 1) && (Q[i] == 0)) {",
    "            alpha_0 = alpha_1_prev + one_zero(ex[i]);",
    "            alpha_1 = alpha_1_prev + one_one(ex[i], e1);",
    "          } else if ((Q[i - 1] == 1) && (Q[i] == 1)) {",
    "            alpha_1 = alpha_1_prev + one_one(ex[i], e1);",
    "          }",
    "        }",
    "      }",
    "    }",
    "  ",
    "    // Return",
    "    real out;",
    "    if (Q[n_year] == 0) {",
    "      out = log_sum_exp(alpha_0, alpha_1);",
    "    } else if (Q[n_year] == 1) {",
    "      out = alpha_1;",
    "    }",
    "    return(out);",
    "  }",
    sep = "\n")
}

##### multi-colex lpmf #####

#' Create Stan code for likelihood function occupancy_multi_colex_lpmf.
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @param max_year Literal integer maximum number of years (or seasons) visited
#'    in any series
#' @return Character string of Stan code corresponding to occupancy_multi_colex_lpmf
#' @noRd
make_occupancy_multi_colex_lpmf <- function (max_rep, max_year) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  assertthat::assert_that(
    is_one_pos_int(max_year, m = 1),
    msg = "max_year must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_multi_colex_lpmf(
    array[] int y, // detection data
    vector mu, // linear predictor for detection
    vector occ, // linear predictor for initial occupancy. Elements after vint1[1] irrelevant.
    vector colo, // linear predictor for colonization. Elements after vint2[1] irrelevant.
    vector ex, // linear predictor for extinction. Elements after vint2[1] irrelevant.
    array[] int vint1, // # of series (# of HMMs). Elements after 1 irrelevant.
    array[] int vint2, // # units (series-years). Elements after 1 irrelevant.
    array[] int vint3, // # years per series. Elements after vint1[1] irrelevant.
    array[] int vint4, // # sampling events per unit (n_rep). Elements after vint2[1] irrelevant.
    array[] int vint5, // Indicator for > 0 detections (Q). Elements after vint2[1] irrelevant.
  
  // indices for jth unit (first rep) for each series. Elements after vint1[1] irrelevant."
  
  sf_text2.1 <- paste0("    array[] int vint", 5 + (1:max_year), collapse = ",\n")
  
  sf_text2.2 <- ",\n"
  
  sf_text2 <- paste0(sf_text2.1, sf_text2.2)
  
  sf_text3 <- "// indices for jth repeated sampling event to each unit (elements after vint2[1] irrelevant):"
  
  sf_text4 <- paste0("    array[] int vint", 5 + max_year + (1:max_rep), collapse = ",\n")
  
  sf_text5 <- paste0(") {
  // Create array of the unit indices that correspond to each series.
    array[vint1[1], ", max_year, "] int unit_index_array;")
  
  sf_text6.1 <- "      unit_index_array[,"
  sf_text6.2 <- 1:max_year
  sf_text6.3 <- "] = vint"
  sf_text6.4 <- 5 + (1:max_year)
  sf_text6.5 <- "[1:vint1[1]];\n"
  sf_text6 <- paste0(sf_text6.1, sf_text6.2, sf_text6.3, sf_text6.4, sf_text6.5, collapse = "")
  
  sf_text7 <- paste0("
  // Create array of the rep indices that correspond to each unit.
    array[vint2[1], ", max_rep, "] int visit_index_array;")
  
  sf_text8.1 <- "      visit_index_array[,"
  sf_text8.2 <- 1:max_rep
  sf_text8.3 <- "] = vint"
  sf_text8.4 <- 5 + max_year + (1:max_rep)
  sf_text8.5 <- "[1:vint2[1]];\n"
  sf_text8 <- paste0(sf_text8.1, sf_text8.2, sf_text8.3, sf_text8.4, sf_text8.5, collapse = "")
  
  sf_text9 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    for (i in 1:vint1[1]) {
      int n_year = vint3[i];
      array[n_year] int Q = vint5[unit_index_array[i,1:n_year]];
      array[n_year] int n_obs = vint4[unit_index_array[i,1:n_year]];
      int max_obs = max(n_obs);
      array[n_year, max_obs] int y_i;
      real occ_i = occ[unit_index_array[i,1]];
      vector[n_year] colo_i = to_vector(colo[unit_index_array[i,1:n_year]]);
      vector[n_year] ex_i = to_vector(ex[unit_index_array[i,1:n_year]]);
      array[n_year] row_vector[max_obs] det_i;
      
      for (j in 1:n_year) {
        if (n_obs[j] > 0) {
          y_i[j, 1:n_obs[j]] = y[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]];
          det_i[j, 1:n_obs[j]] = to_row_vector(mu[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]]);
        }
      }
      lp += forward_colex(n_year, Q, n_obs, y_i, occ_i, colo_i, ex_i, det_i);
    }
    return(lp);
  }
"
  
  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sf_text6, 
               sf_text7, sf_text8, sf_text9, sep = "\n")
  out
}

##### multi_colex_eq_lpmf #####

#' Create Stan code for likelihood function occupancy_multi_colex_eq_lpmf.
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @param max_year Literal integer maximum number of years (or seasons) visited
#'    in any series
#' @return Character string of Stan code corresponding to occupancy_multi_colex_eq_lpmf
#' @noRd
make_occupancy_multi_colex_eq_lpmf <- function (max_rep, max_year) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  assertthat::assert_that(
    is_one_pos_int(max_year, m = 1),
    msg = "max_year must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_multi_colex_eq_lpmf(
    array[] int y, // detection data
    vector mu, // linear predictor for detection
    vector colo, // linear predictor for colonization. Elements after vint2[1] irrelevant.
    vector ex, // linear predictor for extinction. Elements after vint2[1] irrelevant.
    array[] int vint1, // # of series (# of HMMs). Elements after 1 irrelevant.
    array[] int vint2, // # units (series-years). Elements after 1 irrelevant.
    array[] int vint3, // # years per series. Elements after vint1[1] irrelevant.
    array[] int vint4, // # sampling events per unit (n_rep). Elements after vint2[1] irrelevant.
    array[] int vint5, // Indicator for > 0 detections (Q). Elements after vint2[1] irrelevant.
  
  // indices for jth unit (first rep) for each series. Elements after vint1[1] irrelevant."
  
  sf_text2.1 <- paste0("    array[] int vint", 5 + (1:max_year), collapse = ",\n")
  
  sf_text2.2 <- ",\n"
  
  sf_text2 <- paste0(sf_text2.1, sf_text2.2)
  
  sf_text3 <- "// indices for jth repeated sampling event to each unit (elements after vint2[1] irrelevant):"
  
  sf_text4 <- paste0("    array[] int vint", 5 + max_year + (1:max_rep), collapse = ",\n")
  
  sf_text5 <- paste0(") {
  // Create array of the unit indices that correspond to each series.
    array[vint1[1], ", max_year, "] int unit_index_array;")
  
  sf_text6.1 <- "      unit_index_array[,"
  sf_text6.2 <- 1:max_year
  sf_text6.3 <- "] = vint"
  sf_text6.4 <- 5 + (1:max_year)
  sf_text6.5 <- "[1:vint1[1]];\n"
  sf_text6 <- paste0(sf_text6.1, sf_text6.2, sf_text6.3, sf_text6.4, sf_text6.5, collapse = "")
  
  sf_text7 <- paste0("
  // Create array of the rep indices that correspond to each unit.
    array[vint2[1], ", max_rep, "] int visit_index_array;")
  
  sf_text8.1 <- "      visit_index_array[,"
  sf_text8.2 <- 1:max_rep
  sf_text8.3 <- "] = vint"
  sf_text8.4 <- 5 + max_year + (1:max_rep)
  sf_text8.5 <- "[1:vint2[1]];\n"
  sf_text8 <- paste0(sf_text8.1, sf_text8.2, sf_text8.3, sf_text8.4, sf_text8.5, collapse = "")
  
  sf_text9 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    for (i in 1:vint1[1]) {
      int n_year = vint3[i];
      array[n_year] int Q = vint5[unit_index_array[i,1:n_year]];
      array[n_year] int n_obs = vint4[unit_index_array[i,1:n_year]];
      int max_obs = max(n_obs);
      array[n_year, max_obs] int y_i;
      real occ_i = logit(
        inv_logit(colo[unit_index_array[i,1]]) / 
          (inv_logit(colo[unit_index_array[i,1]]) + inv_logit(ex[unit_index_array[i,1]])));
      vector[n_year] colo_i = to_vector(colo[unit_index_array[i,1:n_year]]);
      vector[n_year] ex_i = to_vector(ex[unit_index_array[i,1:n_year]]);
      array[n_year] row_vector[max_obs] det_i;
      
      for (j in 1:n_year) {
        if (n_obs[j] > 0) {
          y_i[j, 1:n_obs[j]] = y[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]];
          det_i[j, 1:n_obs[j]] = to_row_vector(mu[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]]);
        }
      }
      lp += forward_colex(n_year, Q, n_obs, y_i, occ_i, colo_i, ex_i, det_i);
    }
    return(lp);
  }
"
  
  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sf_text6, 
               sf_text7, sf_text8, sf_text9, sep = "\n")
  return(out)
}

##### multi_autologistic_lpmf #####

#' Create Stan code for likelihood function occupancy_multi_autologistic_lpmf.
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @param max_year Literal integer maximum number of years (or seasons) visited
#'    in any series
#' @return Character string of Stan code corresponding to occupancy_multi_autologistic_lpmf
#' @noRd
make_occupancy_multi_autologistic_lpmf <- function (max_rep, max_year) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  assertthat::assert_that(
    is_one_pos_int(max_year, m = 1),
    msg = "max_year must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_multi_autologistic_lpmf(
    array[] int y, // detection data
    vector mu, // linear predictor for detection
    vector occ, // linear predictor for initial occupancy. Elements after vint1[1] irrelevant.
    vector colo, // linear predictor for colonization. Elements after vint2[1] irrelevant.
    vector autologistic, // logit-scale offset for persistence. Elements after vint2[1] irrelevant.
    array[] int vint1, // # of series (# of HMMs). Elements after 1 irrelevant.
    array[] int vint2, // # units (series-years). Elements after 1 irrelevant.
    array[] int vint3, // # years per series. Elements after vint1[1] irrelevant.
    array[] int vint4, // # sampling events per unit (n_rep). Elements after vint2[1] irrelevant.
    array[] int vint5, // Indicator for > 0 detections (Q). Elements after vint2[1] irrelevant.
  
  // indices for jth unit (first rep) for each series. Elements after vint1[1] irrelevant."
  
  sf_text2.1 <- paste0("    array[] int vint", 5 + (1:max_year), collapse = ",\n")
  
  sf_text2.2 <- ",\n"
  
  sf_text2 <- paste0(sf_text2.1, sf_text2.2)
  
  sf_text3 <- "// indices for jth repeated sampling event to each unit (elements after vint2[1] irrelevant):"
  
  sf_text4 <- paste0("    array[] int vint", 5 + max_year + (1:max_rep), collapse = ",\n")
  
  sf_text5 <- paste0(") {
  // Create array of the unit indices that correspond to each series.
    array[vint1[1], ", max_year, "] int unit_index_array;")
  
  sf_text6.1 <- "      unit_index_array[,"
  sf_text6.2 <- 1:max_year
  sf_text6.3 <- "] = vint"
  sf_text6.4 <- 5 + (1:max_year)
  sf_text6.5 <- "[1:vint1[1]];\n"
  sf_text6 <- paste0(sf_text6.1, sf_text6.2, sf_text6.3, sf_text6.4, sf_text6.5, collapse = "")
  
  sf_text7 <- paste0("
  // Create array of the rep indices that correspond to each unit.
    array[vint2[1], ", max_rep, "] int visit_index_array;")
  
  sf_text8.1 <- "      visit_index_array[,"
  sf_text8.2 <- 1:max_rep
  sf_text8.3 <- "] = vint"
  sf_text8.4 <- 5 + max_year + (1:max_rep)
  sf_text8.5 <- "[1:vint2[1]];\n"
  sf_text8 <- paste0(sf_text8.1, sf_text8.2, sf_text8.3, sf_text8.4, sf_text8.5, collapse = "")
  
  sf_text9 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    vector[size(colo)] ex = - (colo + autologistic);
    for (i in 1:vint1[1]) {
      int n_year = vint3[i];
      array[n_year] int Q = vint5[unit_index_array[i,1:n_year]];
      array[n_year] int n_obs = vint4[unit_index_array[i,1:n_year]];
      int max_obs = max(n_obs);
      array[n_year, max_obs] int y_i;
      real occ_i = occ[unit_index_array[i,1]];
      vector[n_year] colo_i = to_vector(colo[unit_index_array[i,1:n_year]]);
      vector[n_year] ex_i = to_vector(ex[unit_index_array[i,1:n_year]]);
      array[n_year] row_vector[max_obs] det_i;
      
      for (j in 1:n_year) {
        if (n_obs[j] > 0) {
          y_i[j, 1:n_obs[j]] = y[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]];
          det_i[j, 1:n_obs[j]] = to_row_vector(mu[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]]);
        }
      }
      lp += forward_colex(n_year, Q, n_obs, y_i, occ_i, colo_i, ex_i, det_i);
    }
    return(lp);
  }
"
  
  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sf_text6, 
               sf_text7, sf_text8, sf_text9, sep = "\n")
  return(out)
}


##### multi_autologistic_eq_lpmf #####

#' Create Stan code for likelihood function occupancy_multi_autologistic_eq_lpmf.
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @param max_year Literal integer maximum number of years (or seasons) visited
#'    in any series
#' @return Character string of Stan code corresponding to occupancy_multi_autologistic_eq_lpmf
#' @noRd
make_occupancy_multi_autologistic_eq_lpmf <- function (max_rep, max_year) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  assertthat::assert_that(
    is_one_pos_int(max_year, m = 1),
    msg = "max_year must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_multi_autologistic_eq_lpmf(
    array[] int y, // detection data
    vector mu, // linear predictor for detection
    vector colo, // linear predictor for colonization. Elements after vint2[1] irrelevant.
    vector autologistic, // logit-scale offset for persistence. Elements after vint2[1] irrelevant.
    array[] int vint1, // # of series (# of HMMs). Elements after 1 irrelevant.
    array[] int vint2, // # units (series-years). Elements after 1 irrelevant.
    array[] int vint3, // # years per series. Elements after vint1[1] irrelevant.
    array[] int vint4, // # sampling events per unit (n_rep). Elements after vint2[1] irrelevant.
    array[] int vint5, // Indicator for > 0 detections (Q). Elements after vint2[1] irrelevant.
  
  // indices for jth unit (first rep) for each series. Elements after vint1[1] irrelevant."
  
  sf_text2.1 <- paste0("    array[] int vint", 5 + (1:max_year), collapse = ",\n")
  
  sf_text2.2 <- ",\n"
  
  sf_text2 <- paste0(sf_text2.1, sf_text2.2)
  
  sf_text3 <- "// indices for jth repeated sampling event to each unit (elements after vint2[1] irrelevant):"
  
  sf_text4 <- paste0("    array[] int vint", 5 + max_year + (1:max_rep), collapse = ",\n")
  
  sf_text5 <- paste0(") {
  // Create array of the unit indices that correspond to each series.
    array[vint1[1], ", max_year, "] int unit_index_array;")
  
  sf_text6.1 <- "      unit_index_array[,"
  sf_text6.2 <- 1:max_year
  sf_text6.3 <- "] = vint"
  sf_text6.4 <- 5 + (1:max_year)
  sf_text6.5 <- "[1:vint1[1]];\n"
  sf_text6 <- paste0(sf_text6.1, sf_text6.2, sf_text6.3, sf_text6.4, sf_text6.5, collapse = "")
  
  sf_text7 <- paste0("
  // Create array of the rep indices that correspond to each unit.
    array[vint2[1], ", max_rep, "] int visit_index_array;")
  
  sf_text8.1 <- "      visit_index_array[,"
  sf_text8.2 <- 1:max_rep
  sf_text8.3 <- "] = vint"
  sf_text8.4 <- 5 + max_year + (1:max_rep)
  sf_text8.5 <- "[1:vint2[1]];\n"
  sf_text8 <- paste0(sf_text8.1, sf_text8.2, sf_text8.3, sf_text8.4, sf_text8.5, collapse = "")
  
  sf_text9 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    vector[size(colo)] ex = - (colo + autologistic);
    for (i in 1:vint1[1]) {
      int n_year = vint3[i];
      array[n_year] int Q = vint5[unit_index_array[i,1:n_year]];
      array[n_year] int n_obs = vint4[unit_index_array[i,1:n_year]];
      int max_obs = max(n_obs);
      array[n_year, max_obs] int y_i;
      real occ_i = logit(
        inv_logit(colo[unit_index_array[i,1]]) / 
          (inv_logit(colo[unit_index_array[i,1]]) + inv_logit(ex[unit_index_array[i,1]])));
      vector[n_year] colo_i = to_vector(colo[unit_index_array[i,1:n_year]]);
      vector[n_year] ex_i = to_vector(ex[unit_index_array[i,1:n_year]]);
      array[n_year] row_vector[max_obs] det_i;
      
      for (j in 1:n_year) {
        if (n_obs[j] > 0) {
          y_i[j, 1:n_obs[j]] = y[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]];
          det_i[j, 1:n_obs[j]] = to_row_vector(mu[visit_index_array[unit_index_array[i, j], 1:n_obs[j]]]);
        }
      }
      lp += forward_colex(n_year, Q, n_obs, y_i, occ_i, colo_i, ex_i, det_i);
    }
    return(lp);
  }
"
  
  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sf_text6, 
               sf_text7, sf_text8, sf_text9, sep = "\n")
  return(out)
}

##### Single threaded #####
# Not currently in use
make_occupancy_single_threaded_lpmf <- function (max_rep, grainsize) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  assertthat::assert_that(
    is_one_pos_int(grainsize),
    msg = "grainsize must be an integer greater than 0"
  )

  
  sf_text0.1 <- "  real occupancy_single_threaded_lpmf(
    array[] int y, // detection data
    vector occ,
    vector mu, // lin pred for detection
    array[] int vint1, // # units (n_unit). Elements after 1 irrelevant.
    array[] int vint2, // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
    array[] int vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
  
  // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):"
  
  sf_text0.2 <- paste0("    array[] int vint", 3 + (1:max_rep), collapse = ",\n")
  
  sf_text0.3 <- ") {"

  sf_text0 <- paste(sf_text0.1, sf_text0.2, sf_text0.3, sep = "\n")
  sf_text1.1 <- 
    "    real lp = reduce_sum(
      // partial sum function
      occupancy_single_partial_sum,
      // occupancy lp to slice
      to_array_1d(occ[1:vint1[1]]),"
  sf_text1.2 <- paste0("      ", grainsize, ", // grainsize")
  sf_text1.3 <- "      y, // detection data
      mu, // lin pred for detection
      vint1, // # units (n_unit). Elements after 1 irrelevant.
      vint2, // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
      vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
  
      // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):"
  
  sf_text1 <- paste(sf_text1.1, sf_text1.2, sf_text1.3, sep = "\n")
  
  sf_text2 <- paste0("      vint", 3 + (1:max_rep), collapse = ",\n")
  
  sf_text3 <- "    );
    return(lp);
  }
"
  
  out <- paste(sf_text0, sf_text1, sf_text2, sf_text3, sep = "\n")
  return(out)
}

#' Create Stan code for likelihood function occupancy_single_partial_sum.
#' @param max_rep Literal integer maximum number of repeated sampling events at 
#'    any unit.
#' @return Character string of Stan code corresponding to occupancy_single_lpmf
#' @noRd
make_occupancy_single_partial_sum <- function (max_rep) {
  assertthat::assert_that(
    is_one_pos_int(max_rep, m = 1),
    msg = "max_rep must be an integer greater than 1"
  )
  
  sf_text1 <- "  real occupancy_single_partial_sum(
    array[] real occ, // sliced linpred for occupancy
    int start,
    int end,
    array[] int y, // detection data
    vector mu, // lin pred for detection
    array[] int vint1, // # units (n_unit). Elements after 1 irrelevant.
    array[] int vint2, // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
    array[] int vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
  
  // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):"
  
  sf_text2 <- paste0("    array[] int vint", 3 + (1:max_rep), collapse = ",\n")
  
  sf_text3 <- paste0(") {
  // Create array of the rep indices that correspond to each unit.
    array[1 + end - start, ", max_rep, "] int index_array;")
  
  sf_text4.1 <- "      index_array[,"
  sf_text4.2 <- 1:max_rep
  sf_text4.3 <- "] = vint"
  sf_text4.4 <- 3 + (1:max_rep)
  sf_text4.5 <- "[start : end];\n"
  sf_text4 <- paste0(sf_text4.1, sf_text4.2, sf_text4.3, sf_text4.4, sf_text4.5, collapse = "")
  
  sf_text5 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    for (i in 1 : (1 + end - start)) {
      array[vint2[i + start - 1]] int indices = index_array[i, 1:vint2[i + start - 1]];
      if (vint3[i + start - 1] == 1) {
        lp += bernoulli_logit_lpmf(1 | occ[i]);
        lp += bernoulli_logit_lpmf(y[indices] | mu[indices]);
      }
      if (vint3[i + start - 1] == 0) {
        lp += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                              sum(log1m_inv_logit(mu[indices])), bernoulli_logit_lpmf(0 | occ[i]));
      }
    }
    return(lp);
  }
"
  
  out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sep = "\n")
  return(out)
}




