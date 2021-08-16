#' Numerically stable log inverse logit
#' @param x logit-scale probability
#' @return the logarithm of the inverse logit of x
log_inv_logit <- function (x) {
  out <- ifelse(x < 0.0, 
                x - log1p(exp(x)), 
                -log1p(exp(-x)))
  return(out)
}


#' Numerically stable log one-minus inverse logit
#' @param x logit-scale probability
#' @return the logarithm of one minus the inverse logit of x
log1m_inv_logit <- function (x) {
  out <- ifelse(x > 0.0, 
                -x - log1p(exp(-x)), 
                -log1p(exp(x)))
  return(out)
}