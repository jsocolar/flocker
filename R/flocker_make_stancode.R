#' Make Stan code for a flocker model with visit-specific detection covariates.
#' Primarily for internal use in \code{flocker()}.
#' @param model_formula Model formula (occupancy + detection)
#' @param flocker_data data, generally the output of \code{make_flocker_data()}.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic effect)
#' @return A character string giving a Stan model.
#' @export
flocker_make_stancode <- function(model_formula, flocker_data, data2){
  data <- flocker_data$flocker_data
  stancode <- brms::make_stancode(model_formula, 
                                  data = data, data2 = data2,
                                  family = brms::bernoulli(link = "logit"))
  # append package info
  stancode <- sub("generated with brms", "generated with flocker 0.1.0 and brms", 
                  stancode)
  # inserts 
  data_insert <- "
  int <lower=1> nvisit[N_occ]; // number of visits to each site
  int <lower=1> max_visit; // maximum number of visits
  int <lower=1> indices[N_occ, max_visit]; // indices for each visit"
  lik_insert <- "
      vector[N_det] mu_det2 = mu_det + Xc_det*b_det;
      for (i in 1:N_occ) {
      int det_indices[nvisit[i]] = indices[i, 1:nvisit[i]];
      if (Y_occ[i] == 1) {
        target += bernoulli_logit_lpmf(1 | mu_occ[i]);
        target += bernoulli_logit_lpmf(Y_det[det_indices] | mu_det2[det_indices]);
      }
      if (Y_occ[i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | mu_occ[i]) + 
          sum(log1m_inv_logit(mu_det2[det_indices])), bernoulli_logit_lpmf(0 | mu_occ[i]));
      }
    }
  }"
  # update stancode 
  stancode1 <- sub("data \\{(.*?)\\}", 
                   paste0("data {\\1", data_insert, "\n}"), 
                   stancode)
  stancode2 <- sub("target += bernoulli_logit_lpmf(Y_occ | mu_occ);", 
                   lik_insert, 
                   stancode1, fixed = TRUE)
  stancode3 <- sub("target \\+\\= bernoulli_logit_glm_lpmf\\(Y_det.*?}", 
                   "", 
                   stancode2)
  stancode3
}
