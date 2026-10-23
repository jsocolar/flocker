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


##### array manipulation #####
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

#' Extract the first column while retaining the draw dimension
#' @param x three-dimensional array with posterior draws in the final dimension
#' @return matrix with first-dimension elements in rows and draws in columns
#' @noRd
first_column_draw_matrix <- function(x) {
  assertthat::assert_that(length(dim(x)) == 3)
  matrix(x[, 1, , drop = FALSE], nrow = dim(x)[1])
}

##### Bookkeeping #####
#' Return the version of the loaded flocker namespace
#' @return character package version
#' @noRd
flocker_version <- function() {
  as.character(unname(getNamespaceVersion("flocker")))
}

#' reconstruct the unit-to-group mapping from packed flocker data
#' @param data packed data frame from a two-level flocker model
#' @return integer vector assigning each unit to a level-two group
#' @noRd
get_unit_group <- function(data) {
  n_group <- data$ff_n_group[1]
  starts <- data$ff_group_start[seq_len(n_group)]
  rep(seq_len(n_group), data$ff_n_unit[starts])
}

#' Reserved names in flocker
#' @return character vector of regexes matching reserved variable names
#' @noRd
flocker_reserved <- function() {
  c("^ff_", "^\\.", "^Intercept$", "^occ$", "^det$", "^colo$", "^ex$", "^autologistic$", "^Omega$")
}

#' Model types in flocker as outputted by `flock`
#' @return character vector of types
#' @noRd
flocker_model_types <- function() {
  c("single", # single-season garden-variety model
    "single_C", # single without rep-varying covariates
    "twolevel_single", # single-season two-level occupancy model
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
    "twolevel_single", # single-season two-level occupancy model
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
    "twolevel_single",
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
      "single", "single_C", "twolevel_single", "augmented",
      "multi", "multi", "multi", "multi"
    ),
    data_input_type = c(
      "single", "single", "twolevel_single", "augmented",
      "multi", "multi", "multi", "multi"
    )
  )
}

##### flocker_fit manipulation #####
#' Test whether an object is of class flocker_fit
#' @param x object to be tested
#' @return logical
#' @noRd
is_flocker_fit <- function(x) {
  inherits(x, "flocker_fit")
}

#' Extract non-data metadata from a flocker data or fit object
#' @param x a flocker_data or flocker_fit object
#' @return a list of metadata
#' @noRd
get_flocker_metadata <- function(x) {
  if (is_flocker_data(x)) {
    if (!is.list(x$metadata)) {
      stop("Incompatible flocker_data structure. Pre-2.0 objects require reformatting with make_flocker_data().")
    }
    return(x$metadata)
  }
  assertthat::assert_that(
    is_flocker_fit(x),
    msg = "x must be a flocker_data or flocker_fit object"
  )
  metadata <- attr(x, "flocker_metadata")
  assertthat::assert_that(
    is.list(metadata),
    msg = paste0(
      "the flocker_fit object does not contain the metadata required for ",
      "this operation"
    )
  )
  metadata
}

# Check the schema even when prediction bypasses get_positions (e.g. a
# plain data frame).
validate_postprocessing_data <- function(flocker_fit, new_data = NULL) {
  validate_packed_schema(flocker_fit)
  if (is_flocker_data(new_data)) {
    validate_packed_schema(new_data)
    if (!identical(get_flocker_metadata(new_data)$type,
                   attr(flocker_fit, "data_type"))) {
      stop("new_data is formatted for a different model type than flocker_fit")
    }
  }
  invisible(NULL)
}

validate_packed_schema <- function(x) {
  if (is_flocker_fit(x)) {
    type <- attr(x, "data_type")
    metadata <- attr(x, "flocker_metadata")
  } else {
    metadata <- get_flocker_metadata(x)
    type <- metadata$type
  }
  if (!identical(metadata$schema_version, 2L)) {
    stop(paste0("Incompatible ", switch(type, single = "single-season", single_C = "rep-constant single-season", multi = "multiseason", twolevel_single = "two-level", augmented = "augmented"),
                " packing schema. Pre-2.0 objects require a compatible 1.x installation or reformatting and refitting."))
  }
  invisible(NULL)
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
  twolevel_single = c("occ", "det", "Omega"),
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
#'   If `TRUE`, returns values associated with each unit in the shape obtained
#'   by selecting the first-visit column of obs. This is relevant in
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
  validate_packed_schema(data_object)
  metadata <- get_flocker_metadata(data_object)
  data_type <- metadata$type
  the_data <- data_object$data

  if (data_type == "single") {
    if (unit_level) {
      return(stats::setNames(metadata$unit_start, metadata$original_dimnames[[1]]))
    }
    return(metadata$original_to_packed)
  }

  if (data_type %in% c("multi", "twolevel_single", "augmented")) {
    return(if (unit_level) metadata$unit_to_packed else metadata$original_to_packed)
  }

  n_unit <- nrow(the_data)
  n_rep <- max(the_data$ff_n_trial)
  index_matrix <- matrix(rep(seq_len(n_unit), n_rep), ncol = n_rep)
  if (!unit_level) index_matrix else index_matrix[, 1]
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
                                  augmented, threads, f_meta = NULL) {
  validate_no_time_series_terms(list(f_occ, f_det, f_col, f_ex, f_auto, f_meta))
  # Check that inputs are valid individually
  validate_params_individually(f_occ, f_det, flocker_data,
                               multiseason, f_col, f_ex, multi_init, f_auto,
                               augmented, threads, f_meta)
  
  metadata <- get_flocker_metadata(flocker_data)

  # Check that parameters are valid in combination
  if (metadata$type == "single") {
    validate_param_combos_single(f_occ, f_det, flocker_data, 
                                 multiseason, f_col, f_ex, multi_init, f_auto,
                                 augmented)
  } else if (metadata$type == "single_C") {
    validate_param_combos_single_C(f_occ, f_det, flocker_data, 
                                   multiseason, f_col, f_ex, multi_init, f_auto,
                                   augmented)
  } else if (metadata$type == "twolevel_single") {
    validate_param_combos_twolevel_single(f_occ, f_det, flocker_data,
                                   multiseason, f_col, f_ex, multi_init, f_auto,
                                   augmented, threads, f_meta)
  } else if (metadata$type == "augmented") {
    validate_param_combos_augmented(f_occ, f_det, flocker_data, 
                                   multiseason, f_col, f_ex, multi_init, f_auto,
                                   augmented, threads, f_meta)
  } else {
    validate_param_combos_multi(f_occ, f_det, flocker_data, 
                                multiseason, f_col, f_ex, multi_init, f_auto,
                                augmented, threads)
  }
  if (brms::is.brmsformula(f_det)) {
    validate_unit_formula_variables(f_det$pforms$occ, f_det$pforms$colo,
                                    f_det$pforms$ex, f_det$pforms$autologistic,
                                    flocker_data)
    validate_meta_formula_variables(f_det$pforms$Omega, flocker_data)
  }
  validate_unit_formula_variables(f_occ, f_col, f_ex, f_auto, flocker_data)
  validate_meta_formula_variables(f_meta, flocker_data)
}

# Reject time-series structures before brms can reorder packed rows. Traverse
# calls rather than text so covariates named e.g. `ar` remain valid. Combined
# brms formulas store subformulas in lists and autocor formulas in attributes.
validate_no_time_series_terms <- function(x) {
  forbidden <- c("ar", "ma", "arma", "cosy", "unstr",
                 "cor_ar", "cor_ma", "cor_arma", "cor_cosy")
  contains_term <- function(x) {
    if (inherits(x, c("cor_arma", "cor_cosy"))) return(TRUE)
    if (is.call(x)) {
      fun <- x[[1L]]
      if (is.call(fun) && as.character(fun[[1L]]) %in% c("::", ":::")) {
        fun <- fun[[3L]]
      }
      if (is.symbol(fun) && as.character(fun) %in% forbidden) return(TRUE)
      if (any(vapply(as.list(x)[-1L], contains_term, logical(1)))) return(TRUE)
    } else if (is.list(x)) {
      if (any(vapply(x, contains_term, logical(1)))) return(TRUE)
    }
    autocor <- attr(x, "autocor", exact = TRUE)
    !is.null(autocor) && contains_term(autocor)
  }
  if (contains_term(x)) {
    stop("Autoregressive and other time-series correlation terms (ar, ma, arma, cosy, unstr) are not supported by flocker.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Check individual validity of params passed to `flock`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_params_individually <- function(f_occ, f_det, flocker_data,
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented, threads, f_meta) {
  # Check that formulas are valid and produce informative errors otherwise
  assertthat::assert_that(
    !brms::is.mvbrmsformula(f_det),
    msg = paste0(
      "mvbrmsformula objects are not supported by flock(). Pass a single ",
      "brmsformula containing all required flocker distributional formulas."
    )
  )
  assertthat::assert_that(
    is_formula(f_det) | brms::is.brmsformula(f_det),
    msg = formula_error("detection")
  )
  assertthat::assert_that(
    is.null(f_occ) | is_formula(f_occ),
    msg = formula_error("occupancy")
  )
  assertthat::assert_that(
    is.null(f_meta) | is_formula(f_meta),
    msg = formula_error("meta-occupancy")
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
    get_flocker_metadata(flocker_data)$type %in% flocker_data_output_types(),
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
  unit_terms <- unique(unlist(lapply(list(f_occ, f_col, f_ex, f_auto), all.vars)))
  assertthat::assert_that(
    !any(unit_terms %in% get_flocker_metadata(flocker_data)$event_covs),
    msg = paste0("Data passed to `make_flocker_data` as event covariates are ",
                 "being used in at least one unit-level formula (i.e. f_occ, ",
                 "f_col, f_ex, f_auto). If these covariates are variable ",
                 "across visits at any unit, this doesn't make sense. If they ",
                 "are constant across visits within every unit, pass these ",
                 "covariates as unit covariates rather than event covariates.")
  )
}

#' Check that meta-occupancy formulas use level-two covariates only
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_meta_formula_variables <- function(f_meta, flocker_data) {
  if (is.null(f_meta)) {
    return(invisible(NULL))
  }
  metadata <- get_flocker_metadata(flocker_data)
  assertthat::assert_that(
    "level2_covs" %in% names(metadata),
    msg = "f_meta is only allowed for two-level models."
  )
  meta_vars <- setdiff(all.vars(f_meta[[length(f_meta)]]), "Intercept")
  assertthat::assert_that(
    !(metadata$level2_group %in% meta_vars),
    msg = paste0(
      "The level-two grouping column (`", metadata$level2_group,
      "`) cannot be used in f_meta because it has only one latent state per level."
    )
  )
  assertthat::assert_that(
    all(meta_vars %in% metadata$level2_covs),
    msg = paste0("All variables in f_meta must be level-two covariates ",
                 "passed in level2_covs.")
  )
}

#' Check validity of some params passed to flock if `type` is `single` or `single_C`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_single_generic <- function(f_occ, f_det, flocker_data, 
                                                 multiseason, f_col, f_ex, multi_init, f_auto,
                                                 augmented) {
  if(!brms::is.brmsformula(f_det)){
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

#' Check validity of params passed to `flock` if `type` is `twolevel_single`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_twolevel_single <- function(f_occ, f_det, flocker_data,
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented, threads, f_meta) {
  if(!brms::is.brmsformula(f_det)){
    assertthat::assert_that(
      is_flocker_formula(f_occ), msg = formula_error("occupancy")
    )
    assertthat::assert_that(
      is_flocker_formula(f_meta), msg = formula_error("meta-occupancy")
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
  assertthat::assert_that(
    all(is.numeric(flocker_data$data$ff_y)) &
      all(flocker_data$data$ff_y %in% c(0, 1)), 
    msg = "All response elements must be 0, 1, or NA"
  )
  assertthat::assert_that(
    is.null(threads),
    msg = "multithreading not supported in two-level single-season models; set threads to NULL"
  )
}

#' Check validity of params passed to `flock` if `type` is `augmented`
#' @inheritParams validate_flock_params
#' @return silent if parameters are valid
#' @noRd
validate_param_combos_augmented <- function(f_occ, f_det, flocker_data, 
                                         multiseason, f_col, f_ex, multi_init, f_auto,
                                         augmented, threads, f_meta) {
  if(!brms::is.brmsformula(f_det)){
    assertthat::assert_that(
      is_flocker_formula(f_occ), msg = formula_error("occupancy")
    )
  }
  assertthat::assert_that(
    is.null(f_meta),
    msg = "f_meta must be NULL for augmented models."
  )
  assertthat::assert_that(
    !brms::is.brmsformula(f_det) || !("Omega" %in% names(f_det$pforms)),
    msg = paste0(
      "Do not include an Omega formula in f_det for augmented models; ",
      "meta-occupancy is internally fixed to ~ 1."
    )
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
  if(!brms::is.brmsformula(f_det)){
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
