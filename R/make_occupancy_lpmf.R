#' Create Stan code for the occupancy_lpmf function
#' Primarily for internal use in \code{flock()}.
#' @param .max_visit Literal integer maximum number of visits to any site.
#' @return Character string of Stan code corresponding to occupancy_lpmf
#' @export


make_occupancy_lpmf <- function (.max_visit) {
  if (!(is.integer(.max_visit) & (.max_visit > 1))) {
    stop(".max_visit must be an integer greater than 1")
  }
  
  sf_text1 <- "real occupancy_lpmf(
  int y[], // detection data
  real mu[], // lin pred for detection
  real occ[], // lin pred for occupancy. Only the first .nsite elements matter.
  int .nsite[], // # sites. Elements after 1 irrelevant.
  int .nvisit[], // # visits per site. Elements after .nsite[1] irrelevant.
  int .Q[], // Indicator for > 0 detections. Elements after .nsite[1] irrelevant.
  
// indices for jth visit to each site (elements after .nsite[1] irrelevant):"
  
  sf_text2 <- paste0("  int .visit_index", 1:.max_visit, "[]", collapse = ",\n")
  
  sf_text3 <- ") {
// Create array of the visit indices that correspond to each site.
  int .index_array[.nsite, .max_visit];
    .index_array[,1] = .visit_index1;
    .index_array[,2] = .visit_index2;"
  
  if (.max_visit > 2) {
    sf_text4.1 <- "    .index_array[,"
    sf_text4.2 <- sf_text4.4 <- 3:.max_visit
    sf_text4.3 <- "] = .visit_index"
    sf_text4.5 <- ";\n"
    sf_text4 <- paste0(sf_text4.1, sf_text4.2, sf_text4.3, sf_text4.4, sf_text4.5, collapse = "")
  } else {
    sf_text4 <- NULL
  }
  
  sf_text5 <- "// Initialize and compute log-likelihood
  .lp <- 0;
  for (i in 1:.nsite[1]) {
    int .indices[.nvisit[i]] = .index_array[i, 1:.nvisit[i]];
    if (Q[i] == 1) {
      .lp += bernoulli_logit_lpmf(1 | occ[i]);
      .lp += bernoulli_logit_lpmf(y[.indices] | mu[.indices]);
    }
    if (Q[i] == 0) {
      .lp += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                            sum(log1m_inv_logit(mu[.indices])), bernoulli_logit_lpmf(0 | occ[i]));
    }
  }
  return(.lp)
}
"

out <- paste(sf_text1, sf_text2, sf_text3, sf_text4, sf_text5, sep = "\n")
return(out)
}
