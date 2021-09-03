#' Numerically stable log inverse logit
#' @param x real number or vector of reals
#' @return the logarithm of the inverse logit of x
#' @export
log_inv_logit <- function (x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  out <- ifelse(x < 0.0, 
                x - log1p(exp(x)), 
                -log1p(exp(-x)))
  return(out)
}

#' Numerically stable log one-minus inverse logit
#' @param x real number or vector of reals
#' @return the logarithm of one minus the inverse logit of x
#' @export
log1m_inv_logit <- function (x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  out <- ifelse(x > 0.0, 
                -x - log1p(exp(-x)), 
                -log1p(exp(x)))
  return(out)
}

#' Test whether an object is of class flocker_fit
#' @param x object to be tested
is_flocker_fit <- function(x) {
  return("flocker_fit" %in% class(x))
}

#' Extract type of flocker_fit from object of class flocker_fit
#' @param x flocker_fit object
type_flocker_fit <- function(x) {
  if (!is_flocker_fit(x)) {
    stop("x must be a flocker_fit object")
  }
  if (!("lik_type" %in% names(attributes(x)))){
    stop("the attributes of x have been altered or corrupted")
  }
  out <- attributes(x)$lik_type
  if (!(out %in% c("C", "V"))) {
    stop("the attributes of x have been altered or corrupted")
  }
  return(out)
}
