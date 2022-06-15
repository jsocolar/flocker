##### Math #####
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

##### array manipulation #####
#' convert matrix to long vector format
#' @param m matrix-like object to expand
#' @return vector
expand_matrix <- function (m) {
  as.vector(as.matrix(m))
}

#' convert array to long matrix format
#' @param a 3D array to expand
#' @return matrix
expand_array_3D <- function(a) {
  out <- a[,,1]
  for (i in 2:dim(a)[3]) {
    out <- rbind(out, a[,,i])
  }
  return(out)
}

#' Get third dimension of an array
#' @param x array for which third dimension is desired
#' @return integer size of third dimension
nslice <- function(x) {
  if (!is.array(x)){
    stop("x is not an array")
  } else if (length(dim(x)) != 3) {
    stop("x is not a 3-D array")
  }
  dim(x)[3]
}

#' rbind a matrix to itself n times
#' @param m a matrix or dataframe
#' @param n number of times to rbind
#' @return matrix
stack_matrix <- function (m, n) {
  do.call(rbind, replicate(n, m, simplify=FALSE))
}

##### Bookkeeping #####
#' Reserved names in flocker
#' @return character vector of regexes matching reserved variable names
flocker_reserved <- function() {
  c("^y$", "^occ_", "^rep_", "^\\.")
}

#' Model types in flocker as outputted by `flock`
#' @return character vector of types
flocker_model_types <- function() {
  c("single", # single-season garden-variety model
    "single_C", # single without rep-varying covariates
    "augmented", # single-season data-augmented multispecies model
    "multi_colex", # multi-season model with explicit colonization/extinction
    "multi_colex_eq", # multi_colex with equilibrium starting probabiltiies
    "multi_autologistic", # multi-season autologistic model
    "single_fp", # single-season false-positive model with known fp probabilities
    "multi_colex_fp", # multi_colex with known fp rates
#    "multi_autologistic_fp", # multi_autologistic with known fp rates
    "multi_colex_eq_fp" # multi_colex_eq with known fp rates
    )
}

#' Options for the \code{type} argument to \code{make_flocker_data}.
#' @return character vector of possible inputs
flocker_data_input_types <- function() {
  c("single",    # covers all model types prefixed with "single"
                 # (rep-constant versus varying is inferred from
                 # existence of event_covs
    "augmented", # the data-augmented single-species model
    "multi"      # covers all model types prefixed with "multi"
    )
}

#' Data types that can be returned by \code{make_flocker_data}.
#' @return character vector of data output types
flocker_data_output_types <- function() {
  c("single",
    "single_C",
    "augmented", 
    "multi"
    )
}

#' flocker_data_output_types that yield rep-constant models
#' @return character vector of data output types
rep_constant_types <- function() {
  c("single_C")
}

#' flocker_data_output_types valid for use with threading
#' @return character vector of data output types
threading_types <- function() {
  c("single_C")
}

#' flocker_data_output_types valid for use with colex models
#' @return character vector of data output types
multi_types <- function() {
  c("multi")
}

#' possible values for the `colex_init` parameter if not NULL
#' @return character vector of options
colex_init_types <- function(){c("explicit", "equilibrium")}

##### flocker_fit manipulation #####
#' Test whether an object is of class flocker_fit
#' @param x object to be tested
#' @return logical
is_flocker_fit <- function(x) {
  return("flocker_fit" %in% class(x))
}

#' Extract lik_type from object of class flocker_fit
#' @param x flocker_fit object
#' @return string giving model type
type_flocker_fit <- function(x) {
  if (!is_flocker_fit(x)) {
    stop("x must be a flocker_fit object")
  }
  if (!("lik_type" %in% names(attributes(x)))){
    stop("the attributes of the flocker_fit object have been altered or corrupted")
  }
  out <- attributes(x)$lik_type
  if (!(out %in% flocker_model_types())) {
    stop("the attributes of the flocker_fit object have been altered or corrupted")
  }
  return(out)
}

#' Get matrix positions corresponding to each row of data in "single" type 
#' flocker_fit
#' @param flocker_fit a "single" type `flocker_fit` object
#' @return an n_row x 2 matrix, where each row contains the indices of the 
#'     corresponding sampling event in the observation dataframe
get_positions_single <- function(flocker_fit) {
  if (attributes(flocker_fit)$lik_type != "single") {
    stop("flocker_fit type is not 'single'")
  }
  n_unit <- flocker_fit$data$n_unit[1]
  index_matrix <- as.matrix(flocker_fit$data[1:n_unit, grepl("^rep_index", 
                                                    names(flocker_fit$data))])
  n_row <- nrow(flocker_fit$data)
  out <- t(sapply(1:n_row, function(x){which(index_matrix == x, arr.ind = T)}))
  return(out)
}

##### Misc #####
#' Check validity of params passed to `flock`
#' @return silent if parameters are valid
validate_flock_params <- function(f_occ, f_det, flocker_data,
                                  multiseason, f_col, f_ex, colex_init, f_auto,
                                  augmented, fp, threads) {
  
  # Check that inputs are valid individually
  validate_params_individually(f_occ, f_det, flocker_data,
                               multiseason, f_col, f_ex, colex_init, f_auto,
                               augmented, fp, threads)
  
  # Check that parameters are valid in combination
  if (flocker_data$type == "single") {
    validate_param_combos_single(f_occ, f_det, flocker_data, 
                                 multiseason, f_col, f_ex, colex_init, f_auto,
                                 augmented, fp)
  } else if (flocker_data$type == "single_C") {
    validate_param_combos_single_C(f_occ, f_det, flocker_data, 
                                   multiseason, f_col, f_ex, colex_init, f_auto,
                                   augmented, fp)
  } else if (flocker_data$type == "augmented") {
    validate_param_combos_augmented(f_occ, f_det, flocker_data, 
                                   multiseason, f_col, f_ex, colex_init, f_auto,
                                   augmented, fp, threads)
  } else {
    assertthat::assert_that(flocker_data$type == "multi") 
    # above line is redundant but included for clarity
    validate_param_combos_multi(f_occ, f_det, flocker_data, 
                                multiseason, f_col, f_ex, colex_init, f_auto,
                                augmented, fp, threads)
  }
}

#' Check individual validity of params passed to `flock`
#' @return silent if parameters are valid
validate_params_individually <- function(f_occ, f_det, flocker_data,
                                         multiseason, f_col, f_ex, colex_init, f_auto,
                                         augmented, fp, threads) {
  # Check that formulas are valid and produce informative errors otherwise
  assertthat::assert_that(
    !is.null(f_det),
    msg = "detection formula must be provided"
  )
  assertthat::assert_that(
    inherits(f_det, "formula"),
    msg = formula_error("detection")
  )
  assertthat::assert_that(
    is.null(f_occ) | inherits(f_occ, "formula"),
    msg = formula_error("occupancy")
  )
  assertthat::assert_that(
    is.null(f_col) | inherits(f_col, "formula"),
    msg = formula_error("colonization")
  )
  assertthat::assert_that(
    is.null(f_ex) | inherits(f_ex, "formula"),
    msg = formula_error("extinction")
  )
  assertthat::assert_that(
    is.null(f_auto) | inherits(f_auto, "formula"),
    msg = formula_error("autologistic")
  )
  
  # Check that flocker_data is valid
  assertthat::assert_that(
    flocker_data$type %in% flocker_data_output_types(),
    msg = paste0("unrecognized data type; flocker_data is corrupt. ",
                 "If this error persists on newly formatted objects produced ",
                 "without errors by `make_flocker_data`, please report a bug.")
  )
  
  # check that multiseason and colex_init are valid
  multiseason2 <- multiseason
  if(is.null(multiseason2)){multiseason2 <- "NULL"}
  colex_init2 <- colex_init
  if(is.null(colex_init2)){colex_init2 <- "NULL"}
  assertthat::assert_that(multiseason2 %in% c("NULL", "colex", "autologistic"))
  assertthat::assert_that(colex_init2 %in% c("NULL", "explicit", "equilibrium"))
  assertthat::assert_that(is.logical(augmented))
  assertthat::assert_that(is.logical(fp))
  
  # Check that threads is valid
  assertthat::assert_that(
    is.null(threads) | 
      ((length(threads) == 1) & (threads > 0) & (threads == as.integer(threads))),
    msg = "threads must be null or a positive integer"
  )
}

#' Check validity of params passed to `flock` if `type` is `single`
#' @return silent if parameters are valid
validate_param_combos_single <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, colex_init, f_auto,
                                         augmented, fp) {
  assertthat::assert_that(
    is_flocker_formula(f_occ), msg = formula_error("occupancy")
  )
  assertthat::assert_that(
    is.null(f_col) & is.null(f_ex) & is.null(f_auto),
    msg = "colonization/extinction/autologistic formulas not allowed in single-season model"
  )
  assertthat::assert_that(
    is.null(multiseason),
    msg = "flocker_data formatted for single season but `multiseason` is not NULL."
  )
  assertthat::assert_that(
    is.null(colex_init),
    msg = "flocker_data formatted for single season but `colex_init` is not NULL."
  )
  assertthat::assert_that(
    isFALSE(augmented),
    msg = paste0("flocker_data not formatted for augmented model, but ",
                 "`augmented` is not FALSE."
    )
  )  
  if(fp) {
    assertthat::assert_that(
      all(is.numeric(flocker_data$data$y)) &
        all(flocker_data$data$y >= 0) &
        all(flocker_data$data$y <= 1),
      msg = paste0("for a single-season fp model, all response elements ",
                   "must be numeric between 0 and 1 inclusive")
    )
  } else {
    assertthat::assert_that(
      all(is.numeric(flocker_data$data$y)) &
        all(flocker_data$data$y %in% c(0, 1)), 
      msg = paste0("for a standard single-season model, all response ",
                   "elements must be 0 or 1."
      )
    )
  }
}

#' Check validity of params passed to `flock` if `type` is `single_C`
#' @return silent if parameters are valid
validate_param_combos_single_C <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, colex_init, f_auto,
                                         augmented, fp) {
  assertthat::assert_that(
    is_flocker_formula(f_occ), msg = formula_error("occupancy")
  )
  assertthat::assert_that(
    is.null(f_col) & is.null(f_ex) & is.null(f_auto),
    msg = "colonization/extinction/autologistic formulas not allowed in single-season model"
  )
  assertthat::assert_that(
    is.null(multiseason),
    msg = "flocker_data formatted for single season but `multiseason` is not NULL."
  )
  assertthat::assert_that(
    is.null(colex_init),
    msg = "flocker_data formatted for single season but `colex_init` is not NULL."
  )
  assertthat::assert_that(
    isFALSE(augmented),
    msg = paste0("flocker_data not formatted for augmented model, but ",
                 "`augmented` is not FALSE."
                )
  )  
  assertthat::assert_that(
    !fp,
    msg = paste0("rep-constant fp likelihood not implemented; reformat data ",
                 "for fp model or set fp to FALSE."
                )
  )
  
  assertthat::assert_that(
    all(flocker_data$data$n_suc >= 0) &
      all(flocker_data$data$n_suc == round(flocker_data$data$n_suc)),
    msg = paste0("flocker_data appers to be corrupt, with `type` set to ",
                 "single_C but response not consisting entirely of ",
                 "nonnegative integers. ",
                 "If this error persists on newly formatted objects produced ",
                 "without errors by `make_flocker_data`, please report a bug."
                 )

  )
}

#' Check validity of params passed to `flock` if `type` is `augmented`
#' @return silent if parameters are valid
validate_param_combos_augmented <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, colex_init, f_auto,
                                         augmented, fp, threads) {
  assertthat::assert_that(
    is_flocker_formula(f_occ), msg = formula_error("occupancy")
  )
  assertthat::assert_that(
    is.null(f_col) & is.null(f_ex) & is.null(f_auto),
    msg = "colonization/extinction/autologistic formulas not allowed in single-season model"
  )
  assertthat::assert_that(
    is.null(multiseason),
    msg = "flocker_data formatted for single season but `multiseason` is not NULL."
  )
  assertthat::assert_that(
    is.null(colex_init),
    msg = "flocker_data formatted for single season but `colex_init` is not NULL."
  )
  assertthat::assert_that(
    isTRUE(augmented),
    msg = paste0("flocker_data formatted for augmented model, but ",
                 "`augmented` is not TRUE."
    )
  )  
  assertthat::assert_that(
    !fp,
    msg = "augmented models not implemented for fp likelihoods."
  )
  assertthat::assert_that(
    all(is.numeric(flocker_data$data$y)) &
      all(flocker_data$data$y %in% c(0, 1)), 
    msg = paste0("for an augmented model, all response ",
                 "elements must be 0 or 1."
    )
  )
  assertthat::assert_that(
    is.null(threads),
    msg("multithreading not supported in augmented models; set threads to NULL")
  )
}

#' Check validity of params passed to `flock` if `type` is `multi`
#' @return silent if parameters are valid
validate_param_combos_multi <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, colex_init, f_auto,
                                         augmented, fp, threads) {
  assertthat::assert_that(
    isFALSE(augmented),
    msg = paste0("flocker_data not formatted for augmented model, but ",
                 "`augmented` is not FALSE."
    )
  ) 
  assertthat::assert_that(multiseason %in% c("colex", "autologistic"))
  if (multiseason == "colex") {
    assertthat::assert_that(
      is_flocker_formula(f_col), msg = formula_error("colonization")
    )
    assertthat::assert_that(
      is_flocker_formula(f_ex), msg = formula_error("extinction")
    )
    assertthat::assert_that(colex_init %in% c("explicit", "equilibrium"))
    if (colex_init == "explicit") {
      assertthat::assert_that(
        is_flocker_formula(f_occ), msg = formula_error("occupancy")
      )
    } else {
      assertthat::assert_that(
        is.null(f_occ), msg = "f_occ must be NULL for equilibrium initialization"
      )
    }
  } else {
    assertthat::assert_that(
      is.null(f_col) & is.null(f_ex),
      msg = "colonization/extinction formulas not allowed in autologistic model"
    )
    assertthat::assert_that(
      is_flocker_formula(f_occ),
      msg = formula_error("occupancy")
    )
    assertthat::assert_that(
      is_flocker_formula(f_auto),
      msg = formula_error("autologistic")
    )
    assertthat::assert_that(
      is.null(colex_init),
      msg = "colex_init must be null in autologistic model"
    )
    assertthat::assert_that(
      !fp, msg = "fp likelihoods not yet implemented in autologistic models"
    )
  }
  y_nonmissing <- flocker_data$data$y[flocker_data$data$y != -99]
  if(fp) {
    assertthat::assert_that(
      all(is.numeric(y_nonmissing)) &
        all(y_nonmissing >= 0) &
        all(y_nonmissing <= 1),
      msg = "all fp probabilities must be between 0 and 1 inclusive"
    )
    if (all(y_nonmissing %in% c(0, 1))) {
      warning(paste0("fp is set to TRUE, but all probabilities are either ",
                     "0 or 1. Setting fp = FALSE yields the same model but ",
                     "should fit faster.")
              )
    }
  } else {
    assertthat::assert_that(
      all(is.numeric(y_nonmissing)) &
        all(y_nonmissing %in% c(0, 1)), 
      msg = "in a multi-season model, all non-missing data must be 0 or 1."
    )
  }
  assertthat::assert_that(
    is.null(threads),
    msg("multithreading not supported in multiseason models; set threads to NULL")
  )
}

#' Create error message for formula with incorrect syntax
#' @param form string giving name of formula
#' @return string giving text of error
formula_error <- function(form) {
  strwrap(paste0("Formula error: ", form, " formula has incorrect syntax."))
}

#' check that object is a formula
#' @param x object to test
#' @return logical; TRUE if x is a formula
is_formula <- function (x) {
  inherits(x, "formula")
}

#' check that object is a valid flocker formula
#' @param x object to test
#' @return logical; TRUE if x is a valid flocker formula
is_flocker_formula <- function (x) {
  is_formula(x) & as.character(x)[1] == "~" & length(x) == 2
}

#' check that object is a named list with no duplicate names
#' @param x object to test
#' @return logical; TRUE if x is a named list with no duplicate names
is_named_list <- function (x) {
  is.list(x) & (length(names(x)) == length(x)) & (!any(duplicated(names(x))))
}

#' get shared elements between two vectors
#' @param x vector
#' @param y vector
#' @return all unique elements that occur at least once in both x and y
shared_elements <- function (x, y) {
  ux <- unique(x)
  ux[ux %in% y]
}

#' get largest index in a vector that is not NA
#' @param x vector to test
#' @param treat_m99_NA logical; should -99 be treated as NA?
#' @return 0 if all elements are NA; otherwise the largest index that is not NA
max_position_not_na <- function (x, treat_m99_NA = FALSE) {
  assertthat::assert_that(
    identical(x, as.vector(x)),
    msg = "x is not a vector"
  )
  if (treat_m99_NA) {
    x[x == -99] <- NA
  }
  if (sum(!is.na(x)) == 0) {
    return(0)
  } else {
    max(which(!is.na(x)))
  }
}
