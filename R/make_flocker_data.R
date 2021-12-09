#' Format data for \code{flocker()}.
#' @param obs An I x J matrix-like object where closure is assumed across rows 
#'  and columns are repeated sampling events. Allowable values are 1 (detection), 
#'  0 (no detection), and NA (no sampling event).
#   The data must be formatted so that all NAs are trailing within their rows.
#' @param unit_covs A dataframe of covariates for each unit that are constant 
#'            across repeated sampling events within closure-units.
#' @param event_covs A named list of I x J matrices, each one corresponding to a covariate
#' that varies across repeated sampling events within closure-units
#' @return A flocker_data list that can be passed as data to \code{flocker()}.
#' @export

make_flocker_data <- function(obs, unit_covs = NULL, event_covs = NULL) {
  if(length(dim(obs)) != 2) {
    stop("obs must have exactly two dimensions.")
  }
  n_unit <- nrow(obs)
  n_rep <- ncol(obs)
  if (n_rep < 2) {
    stop("obs must contain at least two columns.")
  }
  if (any(!(unique(obs) %in% c(0, 1, NA)))) {
    stop("obs contains values other than 0, 1, NA.")
  }
  if (any(is.na(obs[ , 1]))) {
    stop("obs has NAs in its first column.")
  }
  if (n_rep > 2) {
    for (j in 2:(n_rep - 1)) {
      the_nas <- is.na(obs[ , j])
      if (any(the_nas)) {
        the_nas2 <- which(the_nas)
        if (any(!is.na(obs[the_nas2, j+1]))) {
          stop(paste0("Some rows of obs have non-trailing NAs."))
        }
      }
    }
  }
  if (all(is.na(obs[ , n_rep]))) {
    stop("The final column of obs contains only NAs.")
  }
  if (!is.null(unit_covs)) {
    if(nrow(unit_covs) != nrow(obs)) {
      stop("Different numbers of rows found for obs and unit_covs.")
    }
    if (any(is.na(unit_covs))) {
      stop("A constant covariate contains missing values.")
    }
    if ("y" %in% names(unit_covs)) {
      stop(paste0("'y' is a reserved variable name in flocker and is not allowed",
                  "as the name of a constant covariate"))
    }
    if (any(grepl("^\\.", names(unit_covs)))) {
      stop(paste0("variable names starting with '.' are reserved in flocker and ",
                  "are not allowed as covariate names in unit_covs"))
    }
    if (any(grepl("^occ_", names(unit_covs)))) {
      stop(paste0("variable names starting with 'occ_' are reserved in flocker and ",
                  "are not allowed as covariate names in unit_covs"))
    }
    if (any(grepl("^rep_", names(unit_covs)))) {
      stop(paste0("variable names starting with 'rep_' are reserved in flocker and ",
                  "are not allowed as covariate names in unit_covs"))
    }
  }
  if (!is.null(event_covs)) {
    if (!is.list(event_covs)) {
      stop("event_covs must be a list or NULL.")
    } else {
      if (is.null(names(event_covs))) {
        stop("If provided, event_covs must be named.")
      }
      if (!is.null(unit_covs)) {
        if (any(names(event_covs) %in% names(unit_covs))) {
          stop("overlapping names detected between event_covs and unit_covs")
        }
      }
      n_event_covs <- length(event_covs)
      missing_covs <- vector()
      for (ec in 1:n_event_covs) {
        if (!all.equal(dim(event_covs[[ec]]), dim(obs))) {
          stop(paste0("Dimension mismatch found between obs and event_covs[[", ec, "]]."))
        }
        missing_covs <- unique(c(missing_covs, which(is.na(event_covs[[ec]]))))
      }
      if (length(missing_covs) > 0) {
        if (!all(missing_covs %in% which(is.na(obs)))) {
          stop(paste0("An event covariate contains missing values ",
                      "at a position where the response is not missing."))
        }
      }
      if ("y" %in% names(event_covs)) {
        stop(paste0("'y' is a reserved variable name in flocker and is not allowed",
                    "as the name of an event covariate"))
      }
      if (any(grepl("^\\.", names(event_covs)))) {
        stop(paste0("variable names starting with '.' are reserved in flocker and ",
                    "are not allowed as covariate names in event_covs"))
      }
      if (any(grepl("^\\.", names(event_covs)))) {
        stop(paste0("variable names starting with '.' are reserved in flocker and ",
                    "are not allowed as covariate names in event_covs"))
      }
      if (any(grepl("^occ_", names(event_covs)))) {
        stop(paste0("variable names starting with 'occ_' are reserved in flocker and ",
                    "are not allowed as covariate names in event_covs"))
      }
      if (any(grepl("^rep_", names(event_covs)))) {
        stop(paste0("variable names starting with 'rep_' are reserved in flocker and ",
                    "are not allowed as covariate names in event_covs"))
      }
    }
  }
  
  max_rep <- ncol(obs)
  if (is.null(event_covs)) {
    n_trial <- rowSums(!is.na(obs))
    n_suc <- rowSums(obs, na.rm = T)
    flocker_data <- data.frame(n_suc = n_suc, n_trial = n_trial)
    if (!is.null(unit_covs)) {
      flocker_data <- cbind(flocker_data, unit_covs)
    }
    out <- list(data = flocker_data, max_rep = max_rep, 
                type = "C")  # C for unit covariates
  } else {
    flocker_data <- data.frame(y = expand_matrix(obs))
    if (!is.null(unit_covs)) {
      unit_covs_stacked <- 
        do.call(rbind, replicate(max_rep, unit_covs, simplify=FALSE))
      flocker_data <- cbind(flocker_data, unit_covs_stacked)
    }
    if (!is.null(event_covs)) {
      event_covs <- as.data.frame(lapply(event_covs, expand_matrix))
      flocker_data <- cbind(flocker_data, event_covs)
    }
    
    
    flocker_data$n_unit <- c(nrow(obs), 
                             rep(-99, nrow(obs) - 1))
    
    flocker_data$n_rep <- c(matrixStats::rowSums2(!is.na(obs)), 
                            rep(-99, nrow(obs) * (max_rep - 1)))
    
    
    flocker_data$Q <- c(as.integer(matrixStats::rowSums2(obs, na.rm = T) > 0),
                        rep(-99, nrow(obs) * (max_rep - 1)))
    
    # Prepare to add rep indices, and trim flocker_data to existing observations
    flocker_data$unit <- 1:nrow(obs)
    
    rep_indices <- as.data.frame(matrix(data = -99, nrow = nrow(flocker_data),
                                        ncol = max_rep))
    names(rep_indices) <- paste0("rep_index", 1:max_rep)
    
    is_not_na <- !is.na(flocker_data$y)
    rep_index_vec <- rep(-99, max_rep*nrow(obs))
    rep_index_vec[is_not_na] <- cumsum(is_not_na)[is_not_na]
    rep_indices[1:nrow(obs),] <- rep_index_vec
    
    flocker_data <- flocker_data[!is.na(flocker_data$y), ]
    rep_indices <- rep_indices[is_not_na,]
    flocker_data <- cbind(flocker_data, rep_indices)
    
    out <- list(data = flocker_data, max_rep = max_rep,
                type = "V")
  }
  class(out) <- c("list", "flocker_data")
  out
}


# helper function to convert matrix to long vector format
# @param m matrix-like object to expand

expand_matrix <- function (m) {
  out <- as.vector(as.matrix(m))
}
