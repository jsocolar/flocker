% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_flocker_data.R
\name{standard_mfd_checks}
\alias{standard_mfd_checks}
\title{input checking for make_flocker_data}
\usage{
standard_mfd_checks(
  obs,
  unit_covs,
  event_covs,
  type,
  n_aug,
  quiet,
  newdata_checks
)
}
\arguments{
\item{obs}{If \code{type = "single"}, an I x J matrix-like object where 
closure is assumed across rows and columns are repeated sampling events. 
  If \code{type = "multi"}, an I x J x K array where rows are sites or 
species-sites, columns are repeated sampling events, and slices along the 
third dimension are seasons. Allowable values are 1 (detection), 0 (no 
detection), and NA (no sampling event).
   If \code{type = "augmented"}, an L x J x K array where rows L are sites, 
columns J are repeat sampling events, and slices K are species. 
   The data must be packed so that, for a given unit (site, site-species, 
site-timestep, site-species-timestep) all realized visits come before any 
missing visits (NAs are trailing within their rows).}

\item{unit_covs}{If \code{type = "single"} a dataframe of covariates for each 
closure-unit that are constant across repeated sampling events within units.
  If \code{type = "multi"}, a list of such dataframes, one per timestep. All 
dataframes must have identical column names and types, and all
dataframes must have I rows.
  If \code{type = "augmented"}, a dataframe of covariates for each site that
are constant across repeated sampling events within sites (no dependence on
species is allowed).}

\item{event_covs}{If \code{type = "single"}, a named list of I x J matrices, 
each one corresponding to a covariate that varies across repeated sampling 
events within closure-units.
  If \code{type = "multi"}, a named list of I x J x K arrays, each one 
corresponding to a covariate that varies across repeated sampling events 
within closure-units.
  If \code{type = "augmented"}, a named list of L x J matrices, each one
corresponding to a covariate that varies across repeated sampling events
within sites (no dependence on species is allowed).}

\item{type}{The type of occupancy model desired. Options are:
\code{"single"} for a single_season model,
\code{"multi"} for a multi-season (dynamic) model, or
\code{"augmented"} for a single-season multi-species model with 
data-augmentation for never-observed pseudospecies.}

\item{n_aug}{Number of pseudo-species to augment. Only applicable if 
\code{type = "augmented"}.}

\item{quiet}{Hide progress bars and informational messages?}

\item{newdata_checks}{If TRUE, turn off checks that must pass in order
to use the data for model fitting, but not in other contexts (e.g. making
predictions or assessing log-likelihoods over new data).}
}
\description{
input checking for make_flocker_data
}
