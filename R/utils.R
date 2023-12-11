##### Math #####
#' Numerically stable log inverse logit
#' @param x real number or vector of reals
#' @return the logarithm of the inverse logit of x
#' @export
#' @examples
#' log_inv_logit(0)
log_inv_logit <- function (x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  out <- ifelse(x < 0.0, 
                x - log1p(exp(x)), 
                -log1p(exp(-x)))
  out
}

#' Numerically stable log one-minus inverse logit
#' @param x real number or vector of reals
#' @return the logarithm of one minus the inverse logit of x
#' @export
#' @examples
#' log1m_inv_logit(0)
log1m_inv_logit <- function (x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  out <- ifelse(x > 0.0, 
                -x - log1p(exp(-x)), 
                -log1p(exp(x)))
  return(out)
}


#' extended binomial RNG
#' @param p binomial probability
#' @param s number of successes
#' @return a sample representing a possible number of trials
#' @details random samples `y` from the distribution proportional to 
#'   binomial_pmf(y | n, s)
#' @noRd
extended_binomial_rng <- function(p, s) {
  assertthat::assert_that(is.numeric(p) & length(p) == 1)
  assertthat::assert_that(is_one_pos_int(s))
  assertthat::assert_that(p > 0 & p <= 1)
  
  r <- stats::runif(length(p))
  c <- 0
  y <- s - 1
  while(c < r){
    y <- y + 1
    c <- c + p * stats::dbinom(s, y, p)
  }
  y
}

##### array manipulation #####
#' convert matrix-like object to long vector format
#' @param m matrix-like object to expand
#' @return vector
#' @noRd
expand_matrix <- function (m) {
  assertthat::assert_that(
    isTRUE(length(dim(m)) == 2),
    msg = "error in expand_matrix: dimension of input is not 2"
    )
  as.vector(as.matrix(m))
}

#' convert array to long matrix format
#' @param a 3D array to expand
#' @return matrix
#' @noRd
expand_array_3D <- function(a) {
  assertthat::assert_that(
    nslice(a) > 1,
    msg = "error in expand_array_3D: input contains fewer than two slices"
  )
  out <- a[,,1]
  for (i in 2:dim(a)[3]) {
    out <- rbind(out, a[,,i])
  }
  out
}

#' Get third dimension of a 3D array
#' @param x 3D array for which third dimension is desired
#' @return integer size of third dimension
#' @noRd
nslice <- function(x) {
  assertthat::assert_that(
    is.array(x),
    msg = "error in nslice: input is not an array"
  )
  assertthat::assert_that(
    length(dim(x)) == 3,
    msg = "error in nslice: input is not a 3-D array"
  )
  dim(x)[3]
}

#' rbind a matrix or dataframe to itself n times
#' @param m a matrix or dataframe
#' @param n number of times to rbind
#' @return matrix or dataframe
#' @noRd
stack_matrix <- function (m, n) {
  assertthat::assert_that(
    length(dim(m)) == 2, 
    msg = "input not two-dimensional"
  )
  do.call(rbind, replicate(n, m, simplify=FALSE))
}

#' create a new matrix with the dimensions of an old matrix
#' @param m old matrix
#' @param data passed to the data argument of matrix()
#' @param byrow passed to the byrow argument of matrix()
#' @return new matrix
#' @noRd
new_matrix <- function(m, data = NA, byrow = FALSE){
  matrix(data, nrow = nrow(m), ncol = ncol(m), byrow = byrow)
}

#' create a new array with the dimensions of an old array
#' @param m old array
#' @param data passed to the data argument of array()
#' @return new array
#' @noRd
new_array <- function(m, data = NA){
  array(data, dim = dim(m))
}

##### Bookkeeping #####
#' column names created in flocker
#' @param n_rep max number of repeat visits
#' @param n_year max length of hmm series in dynamic models
#' @return character vector of column names created internally by 
#'    make_flocker_data
#' @noRd
flocker_col_names <- function(n_rep = NULL, n_year = NULL) {
  out <- c("ff_y", 
    "ff_n_suc", "ff_n_trial", 
    "ff_Q", "ff_n_unit", "ff_n_rep", "ff_unit",
    "ff_n_series", "ff_n_year", "ff_series", "ff_year", "ff_series_year",
    "ff_n_sp", "ff_species", "ff_superQ")
  if(!is.null(n_rep)) {
    out <- c(out, paste0("ff_rep_index", 1:n_rep))
  }
  if(!is.null(n_year)) {
    out <- c(out, paste0("ff_unit_index", 1:n_year))
  }
  out
}

#' Reserved names in flocker
#' @return character vector of regexes matching reserved variable names
#' @noRd
flocker_reserved <- function() {
  c("^ff_", "^\\.", "^occ$", "^det$", "^colo$", "^ex$", "^autologistic$", "^Omega$")
}

#' Model types in flocker as outputted by `flock`
#' @return character vector of types
#' @noRd
flocker_model_types <- function() {
  c("single", # single-season garden-variety model
    "single_C", # single without rep-varying covariates
    "augmented", # single-season data-augmented multispecies model
    "multi_colex", # multi-season model with explicit colonization/extinction
    "multi_colex_eq", # multi_colex with equilibrium starting probabiltiies
    "multi_autologistic", # multi-season autologistic model
    "multi_autologistic_eq" # multi_autologistic with equilibrium starting probabilities
    )
}

#' Options for the \code{type} argument to \code{make_flocker_data}.
#' @return character vector of possible inputs
#' @noRd
flocker_data_input_types <- function() {
  c("single",    # covers all model types prefixed with "single"
                 # (rep-constant versus varying is inferred from
                 # existence of event_covs)
    "augmented", # the data-augmented multispecies model
    "multi"      # covers all model types prefixed with "multi"
    )
}

#' Data types that can be returned by \code{make_flocker_data}.
#' @return character vector of data output types
#' @noRd
flocker_data_output_types <- function() {
  c("single",
    "single_C",
    "augmented", 
    "multi"
    )
}

#' flocker data type lookup
#' @return dataframe giving lookup table for data input types, data output types
#'   and model types
#' @noRd
fdtl <- function(){
  data.frame(
    model_type = flocker_model_types(),
    data_output_type = c(
      "single", "single_C", "augmented", "multi", "multi", "multi", "multi"
    ),
    data_input_type = c(
      "single", "single", "augmented", "multi", "multi", "multi", "multi"
    )
  )
}

#' flocker_data_output_types that yield rep-constant models
#' @return character vector of data output types
#' @noRd
rep_constant_types <- function() {
  c("single_C")
}

#' flocker_data_output_types valid for use with threading
#' @return character vector of data output types
#' @noRd
threading_types <- function() {
  c("single_C")
}

#' flocker_data_output_types valid for use with colex models
#' @return character vector of data output types
#' @noRd
multi_types <- function() {
  c("multi")
}

#' possible values for the `multi_init` parameter if not NULL
#' @return character vector of options
#' @noRd
multi_init_types <- function(){c("explicit", "equilibrium")}

##### flocker_fit manipulation #####
#' Test whether an object is of class flocker_fit
#' @param x object to be tested
#' @return logical
#' @noRd
is_flocker_fit <- function(x) {
  inherits(x, "flocker_fit")
}

#' Extract lik_type from object of class flocker_fit
#' @param x flocker_fit object
#' @return string giving model type
#' @noRd
type_flocker_fit <- function(x) {
  assertthat::assert_that(
    is_flocker_fit(x),
    msg = "x must be a flocker_fit object"
  )
  assertthat::assert_that(
    all(c("data_type") %in% names(attributes(x))),
    msg = "the attributes of the flocker_fit object have been altered or corrupted"
  )
  a <- attributes(x)
  out <- a$data_type
  if(!is.null(a$multiseason)){
    out <- paste(out, a$multiseason, sep = "_")
    if(a$multi_init == "equilibrium"){
      out <- paste(out, "eq", sep = "_")
    }
  }
  
  assertthat::assert_that(
    out %in% flocker_model_types(),
    msg = "the attributes of the flocker_fit object have been altered or corrupted"
  )
  
  out
}

#' A list of distributional parameters by model type
#' @return A list of distributional parameters by model type
#' @noRd
params_by_type <- list(
  single = c("occ", "det"),
  single_C = c("occ", "det"),
  augmented = c("occ", "det", "Omega"),
  multi_colex = c("occ", "colo", "ex", "det"),
  multi_colex_eq = c("colo", "ex", "det"),
  multi_autologistic = c("occ", "colo", "auto", "det"),
  multi_autologistic_eq = c("colo", "auto", "det")
)


#' Get a matrix/array of row numbers from flocker_data$data in the shape of obs
#' @param data_object a flocker_fit object or a flocker_data object
#' @param unit_level logical. If `FALSE`, returns values associated with each visit
#'   in the shape of obs, with NAs for visits that did not occur.
#'   If `TRUE`, returns values associated with each unit in the shape of 
#'   the slice of obs corresponding to the first visit. This is relevant in
#'   multiseason models, where it is possible to have units (i.e. timesteps) 
#'   that are part of the timeseries and have linear predictors for colonization 
#'   etc, but that received no visits. These units are dropped from the 
#'   formatted data if they are trailing (at the end of the series) but must be 
#'   included otherwise. If `TRUE`, the return will have NAs for the trailing 
#'   units that are dropped, but will have row-numbers representing the rows 
#'   giving the unit covariates for the never-visited units if the timesteps are 
#'   not trailing.
#' @return a matrix/array of the same shape as the observations passed to 
#'   `make_flocker_data` (for an augmented model, in the shape of the augmented
#'   data created after augmentation) where each element gives the row number in 
#'   flocker_data$data where the corresponding element resides. If the data
#'   were formatted for a data-augmented model, the returned array will have
#'   extra slices for all of the augmented species.
#' @noRd
get_positions <- function(data_object, unit_level = FALSE) {
  assertthat::assert_that(
    isTRUE(is_flocker_fit(data_object)) | isTRUE(is_flocker_data(data_object)),
    msg = "the data object must either be a flocker_fit or a flocker_data object"
  )
  the_data <- data_object$data
  if(is_flocker_fit(data_object)) {
    data_type <- attributes(data_object)$data_type
  } else {
    data_type <- data_object$type
  }
  
  if(data_type == "single") {
    n_unit <- the_data$ff_n_unit[1]
    n_rep <- max(the_data$ff_n_rep[seq_len(n_unit)], na.rm = TRUE)
    
    index_matrix <- as.matrix(
      the_data[seq_len(n_unit), grepl("^ff_rep_index", names(the_data))]
      )
    assertthat::assert_that(ncol(index_matrix) == n_rep)
    index_matrix[index_matrix == -99] <- NA
    if(!unit_level) {
      return(index_matrix)
    } else {
      return(index_matrix[, 1])
    }
  } else if(data_type == "single_C") {
    n_rows <- nrow(the_data)
    n_cols <- max(the_data$ff_n_trial)
    index_matrix <- matrix(rep(seq_len(n_rows), n_cols), ncol = n_cols)
    if(!unit_level) {
      return(index_matrix)
    } else {
      return(index_matrix[, 1])
    }
  } else if(data_type == "augmented") {
    n_species <- the_data$ff_n_sp[1]
    n_site <- the_data$ff_n_unit[1] / n_species
    max_visit <- max(the_data$ff_n_rep)
    index_array <- array(dim = c(n_site, max_visit, n_species))
    rep_index_frame <- the_data[paste0("ff_rep_index", seq_len(max_visit))]
    for(r in seq_len(nrow(the_data))){
      rep_pos <- which(rep_index_frame == r, arr.ind = TRUE)
      unit_id <- rep_pos[1]
      visit_id <- rep_pos[2]
      sp_id <- the_data$ff_species[r]
      site_id <- unit_id %% n_site
      if(site_id == 0){
        site_id <- n_site
      }
      index_array[site_id, visit_id, sp_id] <- r
    }
    if(!unit_level) {
      return(index_array)
    } else {
      return(as.matrix(index_array[,1,]))
    }
  } else if(data_type == "multi") {
    n_series <- the_data$ff_n_series[1]
    n_unit <- the_data$ff_n_unit[1]
    n_year <- the_data$ff_n_year[seq_len(n_series)]
    n_visit <- the_data$ff_n_rep[seq_len(n_unit)]
    max_year <- max(n_year)
    max_visit <- max(n_visit)
    
    if(!unit_level){
      index_array <- array(dim = c(n_series, max_visit, max_year))
      rep_index_frame <- the_data[paste0("ff_rep_index", seq_len(max_visit))]
      unit_index_frame <- the_data[paste0("ff_unit_index", seq_len(max_year))]
      for(r in seq_len(nrow(the_data))){
        rep_pos <- which(rep_index_frame == r, arr.ind = TRUE)
        unit_id <- rep_pos[1]
        if(!is.na(unit_id)){
          visit_id <- rep_pos[2]
          unit_pos <- which(unit_index_frame == unit_id, arr.ind = TRUE)
          series_id <- unit_pos[1]
          year_id <- unit_pos[2]
          index_array[series_id, visit_id, year_id] <- r
        }
      }
      return(index_array)
    } else {
      index_slice <- array(dim = c(n_series, max_year))
      unit_index_frame <- the_data[paste0("ff_unit_index", seq_len(max_year))]
      for(r in seq_len(nrow(the_data))){
        unit_pos <- which(unit_index_frame == r, arr.ind = TRUE)
        series_id <- unit_pos[1]
        if(!is.na(series_id)){
          year_id <- unit_pos[2]
          index_slice[series_id, year_id] <- r
        }
      }
      return(index_slice)
    }
  }
}

#' Get emission likelihoods given observations and detection probabilities
#' @param state Compute the emission likelihood for absence (0) or presence (1)
#' @param obs A matrix of observation histories. Rows are units, columns are visits.
#' @param det A matrix of detection probabilities
#' @return a vector of emission likelihoods (one per row)
#' @noRd
emission_likelihood <- function(state, obs, det) {
  assertthat::assert_that(is.numeric(state))
  assertthat::assert_that(
    isTRUE(state == 0) | isTRUE(state == 1), 
    msg = "the state must be zero or one"
  )
  assertthat::assert_that(inherits(obs, "matrix"))
  assertthat::assert_that(inherits(det, "matrix"))
  assertthat::assert_that(all(obs >= 0, na.rm = TRUE) & all(obs <= 1, na.rm = TRUE))
  assertthat::assert_that(all(det >= 0, na.rm = TRUE) & all(det <= 1, na.rm = TRUE))
  assertthat::assert_that(identical(dim(obs), dim(det)))
  
  # make sure there are no obs for which det is NA.
  assertthat::assert_that(all(which(is.na(det)) %in% which(is.na(obs))))

  if(state == 0){
    out <- apply(1 - obs, 1, prod, na.rm = TRUE)
  } else {
    out <- apply((1 - obs) * (1 - det) + (obs * det), 1, prod, na.rm = TRUE)
  }
  out
}


#' get Z given emission likelihoods and state probability
#' @param el0 emission likelihood given state 0
#' @param el1 emission likelihood given state 1
#' @param psi_unconditional occupancy probability not conditioned on observation
#' @return occupancy probability conditioned on observation
#' @noRd
Z_from_emission <- function(el0, el1, psi_unconditional){
  assertthat::assert_that(length(el0) == length(el1))
  assertthat::assert_that(length(el0) == length(psi_unconditional))
  
  all_vals <- c(el0, el1, psi_unconditional)
  assertthat::assert_that(!any(is.na(all_vals)))
  assertthat::assert_that(all(all_vals >= 0))
  assertthat::assert_that(all(all_vals <= 1))
  
  out <- psi_unconditional * el1 / 
    (psi_unconditional * el1 + (1 - psi_unconditional) * el0)
  out
}

##### Misc #####
#' Check validity of params passed to `flock`
#' @param f_occ occupancy formula
#' @param f_det detection formula
#' @param flocker_data flocker_data object
#' @param multiseason multiseason flag
#' @param f_col colonization formula
#' @param f_ex extinction formula
#' @param multi_init multi_init flag
#' @param f_auto autologistic formula
#' @param augmented augmented flag
#' @param threads threads number
#' @return silent if parameters are valid
#' @noRd
validate_flock_params <- function(f_occ, f_det, flocker_data,
                                  multiseason, f_col, f_ex, multi_init, f_auto,
                                  augmented, threads) {
  
  # Check that inputs are valid individually
  validate_params_individually(f_occ, f_det, flocker_data,
                               multiseason, f_col, f_ex, multi_init, f_auto,
                               augmented, threads)
  
  # Check that parameters are valid in combination
  if (flocker_data$type == "single") {
    validate_param_combos_single(f_occ, f_det, flocker_data, 
                                 multiseason, f_col, f_ex, multi_init, f_auto,
                                 augmented)
  } else if (flocker_data$type == "single_C") {
    validate_param_combos_single_C(f_occ, f_det, flocker_data, 
                                   multiseason, f_col, f_ex, multi_init, f_auto,
                                   augmented)
  } else if (flocker_data$type == "augmented") {
    validate_param_combos_augmented(f_occ, f_det, flocker_data, 
                                   multiseason, f_col, f_ex, multi_init, f_auto,
                                   augmented, threads)
  } else {
    assertthat::assert_that(flocker_data$type == "multi") 
    # above line is redundant but included for clarity
    validate_param_combos_multi(f_occ, f_det, flocker_data, 
                                multiseason, f_col, f_ex, multi_init, f_auto,
                                augmented, threads)
  }
  validate_unit_formula_variables(f_occ, f_col, f_ex, f_auto, flocker_data)
}

#' Check individual validity of params passed to `flock`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_params_individually <- function(f_occ, f_det, flocker_data,
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented, threads) {
  # Check that formulas are valid and produce informative errors otherwise
  assertthat::assert_that(
    is_formula(f_det) | brms::is.brmsformula(f_det) | brms::is.mvbrmsformula(f_det),
    msg = formula_error("detection")
  )
  assertthat::assert_that(
    is.null(f_occ) | is_formula(f_occ),
    msg = formula_error("occupancy")
  )
  assertthat::assert_that(
    is.null(f_col) | is_formula(f_col),
    msg = formula_error("colonization")
  )
  assertthat::assert_that(
    is.null(f_ex) | is_formula(f_ex),
    msg = formula_error("extinction")
  )
  assertthat::assert_that(
    is.null(f_auto) | is_formula(f_auto),
    msg = formula_error("autologistic")
  )
  
  # Check that flocker_data is valid
  assertthat::assert_that(
    flocker_data$type %in% flocker_data_output_types(),
    msg = paste0("unrecognized data type; flocker_data is corrupt. ",
                 "If this error persists on newly formatted objects produced ",
                 "without errors by `make_flocker_data`, please report a bug.")
  )
  
  # check that multiseason and multi_init are valid
  multiseason2 <- multiseason
  if(is.null(multiseason2)){multiseason2 <- "NULL"}
  multi_init2 <- multi_init
  if(is.null(multi_init2)){multi_init2 <- "NULL"}
  assertthat::assert_that(multiseason2 %in% c("NULL", "colex", "autologistic"))
  assertthat::assert_that(multi_init2 %in% c("NULL", "explicit", "equilibrium"))
  assertthat::assert_that(is.logical(augmented))

  # Check that threads is valid
  assertthat::assert_that(
    is.null(threads) | is_one_pos_int(threads),
    msg = "threads must be null or a positive integer"
  )
}

#' Check that no event covariates are used in unit formulas
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_unit_formula_variables <- function(f_occ, f_col, f_ex, f_auto, flocker_data) {
  # Check that no event covariates are used in unit formulas
  unit_terms <- character(0)
  if(!is.null(f_occ)){
    unit_terms <- c(unit_terms, labels(stats::terms(f_occ)))
  }
  if(!is.null(f_col)){
    unit_terms <- c(unit_terms, labels(stats::terms(f_col)))
  }
  if(!is.null(f_ex)){
    unit_terms <- c(unit_terms, labels(stats::terms(f_ex)))
  }
  if(!is.null(f_auto)){
    unit_terms <- c(unit_terms, labels(stats::terms(f_auto)))
  }
  assertthat::assert_that(
    !any(unit_terms %in% flocker_data$event_covs),
    msg = paste0("Data passed to `make_flocker_data` as event covariates are ",
                 "being used in at least one unit-level formula (i.e. f_occ, ",
                 "f_col, f_ex, f_auto). If these covariates are variable ",
                 "across visits at any unit, this doesn't make sense. If they ",
                 "are constant across visits within every unit, pass these ",
                 "covariates as unit covariates rather than event covariates.")
  )
}

#' Check validity of some params passed to flock if `type` is `single` or `single_C`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_single_generic <- function(f_occ, f_det, flocker_data, 
                                                 multiseason, f_col, f_ex, multi_init, f_auto,
                                                 augmented) {
  if(!(brms::is.brmsformula(f_det) | brms::is.mvbrmsformula(f_det))){
    assertthat::assert_that(
      is_flocker_formula(f_occ), msg = formula_error("occupancy")
    )
  }
  assertthat::assert_that(
    is.null(f_col) & is.null(f_ex) & is.null(f_auto),
    msg = "colonization/extinction/autologistic formulas not allowed in single-season model"
  )
  assertthat::assert_that(
    is.null(multiseason),
    msg = "flocker_data formatted for single season but `multiseason` is not NULL."
  )
  assertthat::assert_that(
    is.null(multi_init),
    msg = "flocker_data formatted for single season but `multi_init` is not NULL."
  )
  assertthat::assert_that(
    isFALSE(augmented),
    msg = paste0("flocker_data not formatted for augmented model, but ",
                 "`augmented` is not FALSE."
    )
  )  
}

#' Check validity of params passed to `flock` if `type` is `single`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_single <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented) {
  validate_param_combos_single_generic(f_occ, f_det, flocker_data, 
                                       multiseason, f_col, f_ex, multi_init, f_auto,
                                       augmented)

  assertthat::assert_that(
    all(is.numeric(flocker_data$data$ff_y)) &
      all(flocker_data$data$ff_y %in% c(0, 1)), 
    msg = "All response elements must be 0, 1, or NA"
  )
}

#' Check validity of params passed to `flock` if `type` is `single_C`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_single_C <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented) {
  validate_param_combos_single_generic(f_occ, f_det, flocker_data, 
   multiseason, f_col, f_ex, multi_init, f_auto,
   augmented)

  assertthat::assert_that(
    all(flocker_data$data$ff_n_suc >= 0) &
      all(flocker_data$data$ff_n_suc == round(flocker_data$data$ff_n_suc)),
    msg = paste0("flocker_data appers to be corrupt, with `type` set to ",
                 "single_C but response not consisting entirely of ",
                 "nonnegative integers. ",
                 "If this error persists on newly formatted objects produced ",
                 "without errors by `make_flocker_data`, please report a bug."
                 )

  )
}

#' Check validity of params passed to `flock` if `type` is `augmented`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_augmented <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented, threads) {
  if(!(brms::is.brmsformula(f_det) | brms::is.mvbrmsformula(f_det))){
    assertthat::assert_that(
      is_flocker_formula(f_occ), msg = formula_error("occupancy")
    )
  }
  assertthat::assert_that(
    is.null(f_col) & is.null(f_ex) & is.null(f_auto),
    msg = "colonization/extinction/autologistic formulas not allowed in single-season model"
  )
  assertthat::assert_that(
    is.null(multiseason),
    msg = "flocker_data formatted for single season but `multiseason` is not NULL."
  )
  assertthat::assert_that(
    is.null(multi_init),
    msg = "flocker_data formatted for single season but `multi_init` is not NULL."
  )
  assertthat::assert_that(
    isTRUE(augmented),
    msg = paste0("flocker_data formatted for augmented model, but ",
                 "`augmented` is not TRUE."
    )
  )  
  assertthat::assert_that(
    all(is.numeric(flocker_data$data$ff_y)) &
      all(flocker_data$data$ff_y %in% c(0, 1)), 
    msg = "All response elements must be 0, 1, or NA"
  )
  assertthat::assert_that(
    is.null(threads),
    msg = "multithreading not supported in augmented models; set threads to NULL"
  )
}

#' Check validity of params passed to `flock` if `type` is `multi`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_multi <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented, threads) {
  assertthat::assert_that(
    isFALSE(augmented),
    msg = paste0("flocker_data not formatted for augmented model, but ",
                 "`augmented` is not FALSE."
    )
  ) 
  assertthat::assert_that(
    isTRUE(multiseason %in% c("colex", "autologistic")),
    msg = "in a multiseason model, `multiseason` must be either 'colex' or 'autologistic'"
  )
  
  assertthat::assert_that(
    isTRUE(multi_init %in% c("explicit", "equilibrium")),
    msg = "in a multiseason model, `multi_init` must be either 'explicit' or 'equilibrium'"
  )
  if(!(brms::is.brmsformula(f_det) | brms::is.mvbrmsformula(f_det))){
    assertthat::assert_that(
      is_flocker_formula(f_col), msg = formula_error("colonization")
    )
    if (multi_init == "explicit") {
      assertthat::assert_that(
        is_flocker_formula(f_occ), msg = formula_error("occupancy")
      )
    } else {
      assertthat::assert_that(
        is.null(f_occ), msg = "f_occ must be NULL for equilibrium initialization"
      )
    }
    if (multiseason == "colex") {
      assertthat::assert_that(
        is_flocker_formula(f_ex), msg = formula_error("extinction")
      )
      assertthat::assert_that(
        is.null(f_auto), msg = "f_auto must be NULL in colex models"
      )
    } else {
      assertthat::assert_that(
        is_flocker_formula(f_auto), msg = formula_error("autologistic")
      )
      assertthat::assert_that(
        is.null(f_ex), msg = "f_ex must be NULL in autologistic models"
      )
    }
  }
  
  y_nonmissing <- flocker_data$data$ff_y[flocker_data$data$ff_y != -99]

  assertthat::assert_that(
    all(is.numeric(y_nonmissing)) &
      all(y_nonmissing %in% c(0, 1)), 
    msg = "All response elements must be 0, 1, or NA."
  )
  assertthat::assert_that(
    is.null(threads),
    msg = "multithreading not supported in multiseason models; set threads to NULL"
  )
}

#' Create error message for formula with incorrect syntax
#' @param form string giving name of formula
#' @return string giving text of error
#' @noRd
formula_error <- function(form) {
  strwrap(paste0("Formula error: ", form, " formula has incorrect syntax."))
}

#' check that object is a formula
#' @param x object to test
#' @return logical; TRUE if x is a formula
#' @noRd
is_formula <- function (x) {
  inherits(x, "formula")
}

#' check that object is a valid flocker formula
#' @param x object to test
#' @return logical; TRUE if x is a valid flocker formula
#' @noRd
is_flocker_formula <- function (x) {
  is_formula(x) & as.character(x)[1] == "~" & length(x) == 2
}

#' check that object is a flocker data object
#' @param x object to test
#' @return logical; TRUE if x is a flocker data object
#' @noRd
is_flocker_data <- function(x) {
  inherits(x, "flocker_data")
}


#' check that object is a named list with no duplicate names
#' @param x object to test
#' @return logical; TRUE if x is a named list with no duplicate names
#' @noRd
is_named_list <- function (x) {
  names <- names(x)
  is.list(x) & 
    length(x) > 0 & 
    (length(names) == length(x)) & 
    (!any(duplicated(names))) &
    !("" %in% names)
}

#' check if an object is a single logical value
#' @param x object to test
#' @return logical; TRUE if x is a single logical value
#' @noRd
is_one_logical <- function (x) {
  identical(x, TRUE) | identical(x, FALSE)
}

#' check if an object is a single integer (or integerish) > m
#' @param x object to test
#' @param m minimum value for x
#' @return logical; TRUE if x is a single positive integer
#' @noRd
is_one_pos_int <- function(x, m = 0) {
  if(!is.numeric(x)){return(FALSE)}
  if(length(x) != 1){return(FALSE)}
  isTRUE(x == floor(x)) & isTRUE(x > m)
}

#' get shared elements between two vectors
#' @param x vector
#' @param y vector
#' @return all unique elements that occur at least once in both x and y
#' @noRd
shared_elements <- function (x, y) {
  ux <- unique(x)
  ux[ux %in% y]
}

#' get largest index in a vector that is not NA
#' @param x vector to test
#' @param treat_m99_NA logical; should -99 be treated as NA?
#' @details if `treat_m99_NA` is `TRUE`, then character `"-99"` will be treated
#'   as `NA` just as numeric `-99` would be.
#' @return 0 if all elements are NA; otherwise the largest index that is not NA
#' @noRd
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

#' remove rownames
#' @param m object whose rownames to remove
#' @return m with rownames removed
#' @noRd
remove_rownames <- function(m) {
  rownames(m) <- NULL
  m
}

#' rbinom without warning about NAs in the probs
#' @param n single integer number of observations
#' @param size number of trials
#' @param prob probability of success on each trial
#' @noRd
rbinom2 <- function(n, size, prob) {
  if(!any(is.na(prob))) {
    out <- stats::rbinom(n, size, prob)
  } else {
    assertthat::assert_that(length(prob) == n)
    assertthat::assert_that(length(size) %in% c(1,n))
    out <- rep(NA, n)
    incl <- !is.na(prob)
    n2 <- sum(incl)
    if(length(size) > 1){
      size2 <- size[incl]
    } else {
      size2 <- size
    }
    prob2 <- prob[incl]
    out[!is.na(prob)] <- stats::rbinom(n2, size2, prob2)
  }
  out
}
