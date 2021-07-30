#' Make Stan code for a flocker model with visit-specific detection covariates.
#' Primarily for internal use in `flocker()`
#' @param f_occ_use Occupancy formula
#' @param f_det_use Detection formula
#' @param flocker_data data, generally the output of `make_flocker_data()`.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic effect)
#' @return A character string giving a Stan model.
#' @export
flocker_make_stancode <- function(model_formula, flocker_data, data2){
  data <- flocker_data$flocker_data
  stancode <- brms::make_stancode(model_formula, 
                                  data = data, data2 = data2,
                                  family = brms::bernoulli(link = "logit"))
  sc <- unlist(strsplit(stancode, "  int prior_only;  // should the likelihood be ignored?", fixed = T))
  sc1 <- gsub("generated with brms", "generated with flocker 0.1.0 and brms", sc[1])
  sc1 <- paste0(sc1, 
                "  int prior_only;  // should the likelihood be ignored?\n",
                "  int <lower=1> nvisit[N_occ]; // number of visits to each site\n",
                "  int <lower=1> max_visit; // maximum number of visits\n",
                "  int <lower=1> indices[N_occ, max_visit]; // indices for each visit\n")
  stancode1 <- paste0(sc1, sc[2])
  
  sc <- unlist(strsplit(stancode1, "target += bernoulli_logit_lpmf(Y_occ | mu_occ);", fixed = T))
  sc1 <- paste0(sc[1],
                "vector[N_det] mu_det2 = mu_det + Xc_det*b_det;\n",
                "  for (i in 1:N_occ) {\n",
                "    int det_indices[nvisit[i]] = indices[i, 1:nvisit[i]];\n",
                "    if (Y_occ[i] == 1) {\n",
                "      target += bernoulli_logit_lpmf(1 | mu_occ[i]);\n",
                "      target += bernoulli_logit_lpmf(Y_det[det_indices] | mu_det2[det_indices]);\n",
                "    }\n",
                "    if (Y_occ[i] == 0) {\n",
                "      target += log_sum_exp(bernoulli_logit_lpmf(1 | mu_occ[i]) + sum(log1m_inv_logit(mu_det2[det_indices])),\n",
                "                            bernoulli_logit_lpmf(0 | mu_occ[i]));\n",
                "    }\n",
                "  }\n",
                "}\n")
  sc_ <- unlist(strsplit(sc[2], "// priors including constants\n"))
  sc2 <- paste0("  // priors including constants\n",
                sc_[2])
  flocker_stancode <- paste0(sc1, sc2)
  flocker_stancode
}
