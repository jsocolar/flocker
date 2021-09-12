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

#' Get matrix positions corresponding to each row of data in rep-varying 
#' flocker_fit
#' @param flocker_fit a rep-varying `flocker_fit` object
#' @return an n_row x 2 matrix, where each row contains the indices of the 
#'     corresponding sampling event in the observation dataframe
get_positions_V <- function(flocker_fit) {
  if (attributes(flocker_fit)$lik_type != "V") {
    stop("flocker_fit type is not 'V'; i.e. the model is not rep-varying")
  }
  n_unit <- flocker_fit$data$n_unit[1]
  index_matrix <- as.matrix(flocker_fit$data[1:n_unit, grepl("^rep_index", 
                                                    names(flocker_fit$data))])
  n_row <- nrow(flocker_fit$data)
  out <- t(sapply(1:n_row, function(x){which(index_matrix == x, arr.ind = T)}))
  return(out)
}
