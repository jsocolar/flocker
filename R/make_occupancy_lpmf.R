#' Create Stan code for occupancy_vv_lpmf with visit-varying covariates
#' Primarily for internal use in \code{flock()}.
#' @param max_visit Literal integer maximum number of visits to any site.
#' @return Character string of Stan code corresponding to occupancy_vv_lpmf
#' @export


make_occupancy_vv_lpmf <- function (max_visit) {
  if (!(is.integer(max_visit) & (max_visit > 1))) {
    stop("max_visit must be an integer greater than 1")
  }
  
  sf_text1 <- "  real occupancy_vv_lpmf(
    int[] y, // detection data
    vector mu, // lin pred for detection
    vector occ, // lin pred for occupancy. Only the first vint1[1] elements matter.
    int[] vint1, // # sites (nsite). Elements after 1 irrelevant.
    int[] vint2, // # visits per site (nvisit). Elements after vint1[1] irrelevant.
    int[] vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
  
  // indices for jth visit to each site (elements after vint1[1] irrelevant):"
  
  sf_text2 <- paste0("    int[] vint", 3 + (1:max_visit), collapse = ",\n")
  
  sf_text3 <- paste0(") {
  // Create array of the visit indices that correspond to each site.
    int index_array[vint1[1], ", max_visit, "];")
  
  sf_text4.1 <- "      index_array[,"
  sf_text4.2 <- 1:max_visit
  sf_text4.3 <- "] = vint"
  sf_text4.4 <- 3 + (1:max_visit)
  sf_text4.5 <- "[1:vint1[1]];\n"
  sf_text4 <- paste0(sf_text4.1, sf_text4.2, sf_text4.3, sf_text4.4, sf_text4.5, collapse = "")
  
  sf_text5 <- "  // Initialize and compute log-likelihood
    real lp = 0;
    for (i in 1:vint1[1]) {
      int indices[vint2[i]] = index_array[i, 1:vint2[i]];
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
return(out)
}



#' Create Stan code occupancy_vc_lpmf with visit-constant covariates
#' Primarily for internal use in \code{flock()}.
#' @return Character string of Stan code corresponding to occupancy_vc_lpmf

make_occupancy_vc_lpmf <- function () {
  "real occupancy_vc_lpmf(int y, real mu, real occ, int trials) {
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
