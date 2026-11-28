# Fixtures are created only when requested. Data-only tests never compile Stan.
windows_cran <- identical(.Platform$OS.type, "windows") &&
  !identical(Sys.getenv("NOT_CRAN"), "true")

fixture_names <- c("single", "single_C", "twolevel", "aug", "multi_colex_ex",
                   "multi_colex_eq", "multi_auto_ex", "multi_auto_eq")

fixture_data <- local({
  memory <- new.env(parent = emptyenv())
  function(name) {
    name <- match.arg(name, fixture_names)
    if (exists(name, memory, inherits = FALSE)) return(memory[[name]])
    value <- withr::with_seed(1, {
      if (name == "twolevel") {
        obs <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), 6, byrow = TRUE)
        make_flocker_data(obs, data.frame(group = factor(rep(letters[1:3], each = 2))),
                          type = "twolevel_single", level2_group = "group", quiet = TRUE)
      } else {
        if (name %in% c("single", "single_C")) {
          fd <- simulate_flocker_data(n_pt = 20, n_sp = if (name == "single") 8 else 5)
          make_flocker_data(fd$obs, fd$unit_covs,
                            if (name == "single") fd$event_covs else NULL, quiet = TRUE)
        } else if (name == "aug") {
          fd <- simulate_flocker_data(augmented = TRUE, n_sp = 10, n_pt = 20)
          dimnames(fd$obs)[[3]] <- paste0("observed_species", seq_len(dim(fd$obs)[3]))
          make_flocker_data(fd$obs, fd$unit_covs, fd$event_covs,
                            type = "augmented", n_aug = 10, quiet = TRUE)
        } else {
          fd <- simulate_flocker_data(n_pt = 20, n_sp = 1, n_season = 6,
            multiseason = if (grepl("colex", name)) "colex" else "autologistic",
            multi_init = if (endsWith(name, "ex")) "explicit" else "equilibrium",
            ragged_rep = TRUE, missing_seasons = TRUE)
          make_flocker_data(fd$obs, fd$unit_covs, fd$event_covs, type = "multi", quiet = TRUE)
        }
      }
    })
    memory[[name]] <- value
    value
  }
})

fixture_fit_args <- function(name) {
  args <- list(flocker_data = fixture_data(name), cores = 2, iter = 8,
               save_warmup = FALSE, refresh = 0, silent = 2, seed = 1,
               backend = "rstan")
  formulas <- switch(name,
    single = list(f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
                  f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | species),
                  chains = 2, warmup = 6),
    single_C = list(f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
                    f_det = ~ 0 + Intercept + uc1 + (1 + uc1 | species)),
    twolevel = list(f_occ = ~ 1, f_det = ~ 1, f_meta = ~ 1, chains = 2, warmup = 6),
    aug = list(f_occ = ~ 0 + Intercept + (1 | ff_species),
               f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
               augmented = TRUE),
    {
      multi <- list(f_det = ~ 0 + Intercept + uc1 + ec1, f_col = ~ 0 + Intercept + uc1,
        multiseason = if (grepl("colex", name)) "colex" else "autologistic",
        multi_init = if (endsWith(name, "ex")) "explicit" else "equilibrium")
      if (endsWith(name, "ex")) multi$f_occ <- ~ 0 + Intercept + uc1
      if (grepl("colex", name)) multi$f_ex <- ~ 0 + Intercept + uc1 else
        multi$f_auto <- ~ 0 + Intercept + uc1
      multi
    })
  # Do not let user brms options change the number of draws in these fixtures.
  if (is.null(formulas$chains)) formulas$chains <- 4
  if (is.null(formulas$warmup)) formulas$warmup <- 4
  # Keep serialized formulas independent of this test environment and its cache.
  formulas <- lapply(formulas, function(x) {
    if (inherits(x, "formula")) environment(x) <- baseenv()
    x
  })
  c(args, formulas)
}

fixture_hash <- function(x) {
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(x, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

fixture_source_signature <- local({
  value <- NULL
  function() {
    if (!is.null(value)) return(value)
    ns <- asNamespace("flocker")
    names <- ls(ns, all.names = TRUE)
    names <- names[!vapply(names, bindingIsActive, logical(1), env = ns)]
    code <- lapply(names, function(name) {
      x <- get(name, ns, inherits = FALSE)
      if (is.function(x)) list(formals = deparse(formals(x)), body = deparse(body(x)))
    })
    stan <- lapply(c("single", "single_C", "multi", "twolevel_single", "augmented"),
      function(name) readLines(system.file("stan", paste0("occupancy_", name, ".stan"),
                                            package = "flocker", mustWork = TRUE)))
    value <<- fixture_hash(list(code = code, stan = stan))
    value
  }
})

fixture_signature <- function(name) {
  args <- fixture_fit_args(name)
  # Formula environments contain test state; only their syntax affects the fit.
  args <- lapply(args, function(x) if (inherits(x, "formula")) deparse(x) else x)
  list(version = 1L, family = name, args = fixture_hash(args),
       source = fixture_source_signature(), R = R.version.string,
       platform = R.version$platform, os = Sys.info()[["sysname"]],
       packages = vapply(c("flocker", "brms", "rstan", "StanHeaders", "Rcpp", "RcppEigen"),
                         function(x) as.character(utils::packageVersion(x)), character(1)))
}

# Cache construction itself is inert. Disk reads, directory creation, and model
# fitting occur only on a cache miss for the particular requested fixture.
new_fixture_cache <- function(build, signature, directory = "") {
  memory <- new.env(parent = emptyenv())
  function(name) {
    if (exists(name, memory, inherits = FALSE)) return(memory[[name]])
    sig <- signature(name)
    path <- if (nzchar(directory)) file.path(directory, paste0("flocker-", name, ".rds")) else NULL
    record <- if (!is.null(path) && file.exists(path)) {
      tryCatch(readRDS(path), error = function(e) NULL)
    } else NULL
    if (is.list(record) && identical(record$signature, sig) && !is.null(record$value)) {
      value <- record$value
    } else {
      value <- build(name)
      if (!is.null(path)) {
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
        tmp <- tempfile(tmpdir = directory, fileext = ".rds")
        on.exit(unlink(tmp))
        saveRDS(list(signature = sig, value = value), tmp)
        if (file.exists(path)) unlink(path)
        if (!file.rename(tmp, path)) stop("Unable to write fixture cache: ", path)
      }
    }
    memory[[name]] <- value
    value
  }
}

fixture_fit <- local({
  get <- new_fixture_cache(
    build = function(name) suppressWarnings(withr::with_seed(1,
      do.call(flock, fixture_fit_args(name)))),
    signature = fixture_signature,
    directory = Sys.getenv("FLOCKER_TEST_CACHE_DIR", unset = ""))
  function(name) {
    name <- match.arg(name, fixture_names)
    if (windows_cran) testthat::skip("Runtime Stan compilation is skipped on CRAN Windows")
    if (!name %in% c("single", "single_C")) testthat::skip_on_cran()
    get(name)
  }
})
