#' Make Stan data for a flocker model with visit-specific detection covariates.
#' Primarily for internal use in \code{flocker()}.
#' @param model_formula The model formula (occupancy + detection)
#' @param flocker_data data, generally the output of \code{make_flocker_data()}.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic effect)
#' @return A standata list.
#' @export
flocker_make_standata <- function(model_formula, flocker_data, data2) {
  data <- flocker_data$flocker_data
  flocker_standata <- brms::make_standata(model_formula, 
                                          data = data, data2 = data2,
                                          family = brms::bernoulli(link = "logit"))
  flocker_standata$nvisit <- flocker_data$.nvisit
  flocker_standata$max_visit <- ncol(flocker_data$.indices)
  flocker_standata$indices <- flocker_data$.indices
  flocker_standata
}
