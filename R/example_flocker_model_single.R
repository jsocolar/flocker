# -------------------------------------------------------------------------
# Why this file exists
#
# The package ships a pre-fitted example model (`example_flocker_model_single`)
# so users can inspect a realistic `brmsfit` object without compiling a Stan
# model themselves. However, serialized `brmsfit` objects can contain internal
# references to namespaces that were loaded when the model was created
# (including packages not listed in flocker's Imports). If such an object is
# stored directly in /data and lazy-loaded, R may attempt to resolve those
# namespaces at install or load time, producing warnings or errors on systems
# where those packages are not installed.
#
# To avoid that, the model is stored as an .rds file in inst/extdata and
# exposed via an active binding. The object is only read from disk on first
# access in an interactive session. During R CMD check (e.g., CRAN/CI), the
# binding instead returns a lightweight placeholder so that package load does
# not depend on optional namespaces embedded in the serialized fit.
# -------------------------------------------------------------------------

.notify_example_model_load <- function() {
  packageStartupMessage(paste0(
    "You are loading example_flocker_model_single, an object supplied by flocker which contains a brmsfit.\n",
    "This brmsfit may contain internal references to packages not required by flocker, including colorspace, lubridate, and terra.\n",
    "If you see an error about these (or other) packages not being available, install the relevant packages to use example_flocker_model_single."
  ))
}

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  
  if (!exists("example_flocker_model_single", envir = ns, inherits = FALSE)) {
    makeActiveBinding(
      "example_flocker_model_single",
      local({
        loaded <- FALSE
        value <- NULL
        warned <- FALSE
        
        function(x) {
          # Setter: forbid assignment
          if (!missing(x)) {
            stop("example_flocker_model_single is read-only.", call. = FALSE)
          }
          
          # During R CMD check, do NOT load the big brmsfit
          if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
            return(flocker_example_unavailable())
          }
          
          # Normal user path: load once on first access
          if (!loaded) {
            
            # Suppress the user-facing load message when the package is being loaded
            # via devtools/pkgload/roxygen (e.g., during devtools::document() or the
            # "Documenting" phase of devtools::check()). In those contexts, pkgload
            # may touch the active binding while building documentation, which would
            # otherwise trigger the message spuriously. We detect this by scanning the
            # current call stack for pkgload/devtools/roxygen2 and only emit the
            # message during normal interactive user access.
            in_pkgload_or_roxygen <- any(grepl(
              "pkgload|devtools|roxygen2",
              vapply(sys.calls(), function(cl) paste(deparse(cl), collapse = ""), character(1))
            ))
            
            if (interactive() && !warned && !in_pkgload_or_roxygen) {
              .notify_example_model_load()
              warned <- TRUE
            }
            
            path <- system.file("extdata", "example_flocker_model_single.rds", package = pkgname)
            if (path == "") {
              stop("Internal example file not found: example_flocker_model_single.rds", call. = FALSE)
            }
            
            value <<- readRDS(path)
            loaded <<- TRUE
          }
          
          value
        }
      }),
      env = ns
    )
  }
}


#' Example fitted flocker model
#'
#' A pre-fitted single-season occupancy model produced by \code{\link{flock}}
#' using simulated data.
#'
#' This object is loaded lazily on first access.
#'
#' @details
#' The object is a \code{brmsfit}. Depending on how it was created, it may
#' contain internal references to optional packages such as
#' \code{colorspace}, \code{lubridate}, or \code{terra}.
#'
#' If you encounter an error about missing namespaces when accessing this
#' object, install the relevant packages.
#'
#' @examples
#' \dontrun{
#' library(flocker)
#' summary(example_flocker_model_single)
#' }
#'
#' @name example_flocker_model_single
#' @export
NULL



flocker_example_unavailable <- function() {
  structure(
    list(message = paste0(
      "example_flocker_model_single is not available during R CMD check.\n",
      "It is a large brmsfit stored in inst/extdata and is intentionally not loaded in checks."
    )),
    class = "flocker_example_unavailable"
  )
}
