#' Format data for \code{flocker()}.
#' @param obs An I x J matrix-like object where closure is assumed across rows (i.e. 
#' rows are sampling sites or species-sites) and columns are visits. Allowable 
#' values are 1 (detection), 0 (no detection), and NA (no visit).
#  The data must be formatted so that all NAs are trailing within their rows.
#' @param constant_covs A dataframe of covariates for each unit that are constant 
#'            across visits.
#' @param visit_covs A named list of I x J matrices, each one corresponding to a covariate
#' that varies by visit.
#' @return A flocker_data list that can be passed as data to \code{flocker()}.
#' @export
make_flocker_data <- function(obs, constant_covs = NULL, visit_covs = NULL) {
  if(length(dim(obs)) != 2) {
    stop("obs must have exactly two dimensions.")
  }
  n_unit <- nrow(obs)
  n_visit <- ncol(obs)
  if (n_visit < 2) {
    stop("obs must contain at least two columns.")
  }
  if (any(!(unique(obs) %in% c(0, 1, NA)))) {
    stop("obs contains values other than 0, 1, NA.")
  }
  if (any(is.na(obs[ , 1]))) {
    stop("obs has NAs in its first column.")
  }
  if (n_visit > 2) {
    for (j in 2:(n_visit - 1)) {
      the_nas <- is.na(obs[ , j])
      if (any(the_nas)) {
        the_nas2 <- which(the_nas)
        if (any(!is.na(obs[the_nas2, j+1]))) {
          stop(paste0("Some rows of obs have non-trailing NAs."))
        }
      }
    }
  }
  if (all(is.na(obs[ , n_visit]))) {
    stop("The final column of obs contains only NAs.")
  }
  if (!is.null(constant_covs)) {
    if(nrow(constant_covs) != nrow(obs)) {
      stop("Different numbers of rows found for obs and constant_covs.")
    }
    if (any(is.na(constant_covs))) {
      stop("A constant covariate contains missing values.")
    }
    if ("y" %in% names(visit_covs)) {
      stop(paste0("'y' is a reserved variable name in flocker and is not allowed",
                  "as the name of a constant covariate"))
    }
    if (any(grepl("^\\.", names(visit_covs)))) {
      stop(paste0("variable names starting with '.' are reserved in flocker and ",
                  "are not allowed as covariate names in constant_covs"))
    }
  }
  if (!is.null(visit_covs)) {
    if (!is.list(visit_covs)) {
      stop("visit_covs must be a list or NULL.")
    } else {
      if (is.null(names(visit_covs))) {
        stop("If provided, visit_covs must be named.")
      }
      if (!is.null(constant_covs)) {
        if (any(names(visit_covs) %in% names(constant_covs))) {
          stop("overlapping names detected between visit_covs and constant_covs")
        }
      }
      n_visit_covs <- length(visit_covs)
      missing_covs <- vector()
      for (vc in 1:n_visit_covs) {
        if (!all.equal(dim(visit_covs[[vc]]), dim(obs))) {
          stop(paste0("Dimension mismatch found between obs and visit_covs[[", vc, "]]."))
        }
        missing_covs <- unique(c(missing_covs, which(is.na(visit_covs[[vc]]))))
      }
      if (length(missing_covs) > 0) {
        if (!all(missing_covs %in% which(is.na(obs)))) {
          stop(paste0("A visit covariate contains missing values ",
                      "at a position where the response is not missing."))
        }
      }
      if ("y" %in% names(visit_covs)) {
        stop(paste0("'y' is a reserved variable name in flocker and is not allowed",
                    "as the name of a visit covariate"))
      }
      
      if (any(grepl("^\\.", names(visit_covs)))) {
        stop(paste0("variable names starting with '.' are reserved in flocker and ",
                    "are not allowed as covariate names in visit_covs"))
      }
    }
  }
  
  max_visit <- ncol(obs)
  if (is.null(visit_covs)) {
    n_trial <- rowSums(!is.na(obs))
    n_suc <- rowSums(obs, na.rm = T)
    flocker_data <- data.frame(n_suc = n_suc, n_trial = n_trial)
    if (!is.null(constant_covs)) {
      flocker_data <- cbind(flocker_data, constant_covs)
    }
    out <- list(data = flocker_data, max_visit = max_visit, 
                type = "C")  # C for visit-constant covariates
  } else {
    flocker_data <- data.frame(y = expand_matrix(obs))
    if (!is.null(constant_covs)) {
      constant_covs_stacked <- 
        do.call(rbind, replicate(max_visit, constant_covs, simplify=FALSE))
      flocker_data <- cbind(flocker_data, constant_covs_stacked)
    }
    if (!is.null(visit_covs)) {
      visit_covs <- as.data.frame(lapply(visit_covs, expand_matrix))
      flocker_data <- cbind(flocker_data, visit_covs)
    }
    
    flocker_data$n_unit <- c(nrow(obs), 
                             rep(-99, nrow(obs) - 1))
    flocker_data$n_visit <- c(apply(obs, 1, function(x){sum(!is.na(x))}), 
                              rep(-99, nrow(obs) * (max_visit - 1)))
    flocker_data$Q <- c(as.integer(rowSums(obs, na.rm = T) > 0),
                         rep(-99, nrow(obs) * (max_visit - 1)))
    
    # Prepare to add visit indices, and trim flocker_data to existing observations
    flocker_data$unit <- 1:nrow(obs)
    flocker_data <- flocker_data[!is.na(flocker_data$y), ]
    visit_indices <- as.data.frame(matrix(data = -99, nrow = nrow(flocker_data),
                                      ncol = max_visit))
    names(visit_indices) <- paste0("visit_index", 1:max_visit)
    for (i in 1:nrow(obs)) {
      visit_indices[i, 1:flocker_data$n_visit[i]] <- which(flocker_data$unit == i)
    }
    flocker_data <- cbind(flocker_data, visit_indices)
    
    out <- list(data = flocker_data, max_visit = max_visit,
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
