#' Format data for `flocker()`.
#' @param obs An I x J matrix-like object where closure is assumed across rows (i.e. 
#' rows are sampling sites or species-sites) and columns are visits. Allowable 
#' values are 1 (detection), 0 (no detection), and NA (no visit).
#  The data must be formatted so that all NAs are trailing within their rows.
#' @param constant_covs A dataframe of covariates for each site (or species-site) that are
#' constant across visits.
#' @param visit_covs A named list of I x J matrices, each one corresponding to a covariate
#' that varies by visit.
#' @return A flocker_data list that can be passed as data to `flocker()`
#' @export
make_flocker_data <- function(obs, constant_covs = NULL, visit_covs = NULL) {
  if(length(dim(obs)) != 2) {
    stop("obs must have exactly two dimensions")
  }
  n_site <- nrow(obs)
  n_visit <- ncol(obs)
  if (n_visit < 2) {
    stop("obs must contain at least two columns")
  }
  if (any(!(unique(obs) %in% c(0, 1, NA)))) {
    stop("obs contains values other than 0, 1, NA")
  }
  if (any(is.na(obs[ , 1]))) {
    stop("obs has NAs in its first column")
  }
  if (n_visit > 2) {
    for (j in 2:(n_visit - 1)) {
      the_nas <- is.na(obs[ , j])
      if (any(the_nas)) {
        the_nas2 <- which(the_nas)
        if (any(!is.na(obs[the_nas2, j+1]))) {
          stop(paste0("some rows of obs have non-trailing NAs"))
        }
      }
    }
  }
  if (all(is.na(obs[ , n_visit]))) {
    warning("the final column of obs contains only NAs")
  }
  if (!is.null(constant_covs)) {
    if(nrow(constant_covs) != nrow(obs)) {
      stop("obs and constant_covs have differing numbers of rows")
    }
  }
  if (!is.null(visit_covs)) {
    if(!is.list(visit_covs)) {
      stop("visit_covs must be provided as a list or must remain NULL")
    } else {
      n_visit_covs <- length(visit_covs)
      for (vc in 1:n_visit_covs) {
        if (!all.equal(dim(visit_covs[[vc]]), dim(obs))) {
          stop(paste0("dimension mismatch found between obs and visit_covs[[", vc, "]]."))
        }
      }
    }
  }
  n_suc <- rowSums(obs, na.rm = T)
  n_trial <- rowSums(!is.na(obs))
  if (is.null(visit_covs)) {
    flocker_data <- data.frame(n_suc = n_suc, n_trial = n_trial)
    if (!is.null(constant_covs)) {
      flocker_data <- cbind(flocker_data, constant_covs)
    }
    out <- list(flocker_data = flocker_data, type = "N")  # N for no visit covariates
  } else {
    occ <- as.numeric(n_suc > 0)
    data1 <- data1_1 <- cbind(occ, as.data.frame(constant_covs))
    for(i in 2:n_visit) {
      data1_1 <- rbind(data1_1, data1)
    }
    data2 <- as.data.frame(c(list(det = as.vector(as.matrix(obs)), lapply(visit_covs, function(x){as.vector(as.matrix(x))}))))
    data2 <- cbind(data2, data1_1)
    data2$.site <- rep(c(1:n_site), n_visit)
    data2$.visit <- rep(c(1:n_visit), each = n_site)
    data2$occ[(n_site + 1):(nrow(data2))] <- -99
    data2 <- data2[!is.na(data2$det), ]
    
    data2$occupancy_subset <- 0
    data2$occupancy_subset[1:n_site] <- 1
    
    
    .indices <- matrix(data = -99, nrow = n_site, ncol = n_visit)
    for (i in 1:nsite) {
      .indices[i, 1:n_trial[i]] <- which(data2$.site == i)
    }
    
    out <- list(flocker_data = data2, .nsite = n_site, .nvisit = n_trial, .indices = .indices, .type = "V")
  }
  class(out <- c("list", "flocker_data"))
  out
}