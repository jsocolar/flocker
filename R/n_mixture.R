# Data-only coefficients for exact N-mixture likelihoods with binomial detection.
# Shared by Poisson, negative-binomial, and zero-inflated Poisson abundance:
# the falling-factorial polynomial depends only on counts, not the mixing law.
# Expand prod(k + j) in falling factorials using
# c'_i = c_(i-1) + (i + j) c_i. Store logs to avoid overflow.
n_mixture_binomial_log_coefficients <- function(y) {
  if (!is.numeric(y) || !length(y) || anyNA(y) ||
      any(!is.finite(y) | y < 0 | y != floor(y)) ||
      sum(y) > .Machine$integer.max - 1) {
    stop("Counts must be a nonempty vector of finite nonnegative integers with a representable total.")
  }
  M <- max(y)
  logc <- 0
  for (count in y[-which.max(y)]) {
    for (j in M - seq_len(count) + 1) {
      degree <- seq_along(logc) - 1
      a <- c(logc + log(degree + j), -Inf)
      b <- c(-Inf, logc)
      hi <- pmax(a, b)
      logc <- hi + log1p(exp(pmin(a, b) - hi))
    }
  }
  logc
}

# Shared count summaries and polynomial storage for those binomial-detection
# families. Offsets/maxima are more general, but this complete helper assumes
# the binomial polynomial expansion; beta-binomial detection cannot reuse it.
# Input is a list of observed count vectors, in packed unit order. Missing
# visits must already have been removed; empty units are not supported here.
# Offsets include a terminal sentinel, permitting slices for singleton units.
n_mixture_binomial_auxiliary_data <- function(counts) {
  if (!is.list(counts) || !length(counts)) stop("Supply a nonempty list of count vectors.")
  coefficients <- lapply(counts, n_mixture_binomial_log_coefficients)
  list(
    ff_nunit = as.integer(length(counts)),
    ff_ncoef = as.integer(sum(lengths(coefficients))),
    ff_unit_offsets = as.integer(c(1, 1 + cumsum(lengths(counts)))),
    ff_coef_offsets = as.integer(c(1, 1 + cumsum(lengths(coefficients)))),
    ff_max_count = as.integer(vapply(counts, max, numeric(1))),
    ff_log_factorial = vapply(counts, function(y) sum(lgamma(y + 1)), numeric(1)),
    ff_log_coef = unlist(coefficients, use.names = FALSE),
    ff_coef_power = unlist(lapply(coefficients, function(x) seq_along(x) - 1),
                           use.names = FALSE)
  )
}

n_mixture_poisson_binomial <- function() {
  brms::custom_family("n_mixture_poisson_binomial", dpars = c("mu", "abund"),
    links = c("identity", "identity"), type = "int", loop = FALSE,
    vars = c("ff_unit_offsets", "ff_coef_offsets", "ff_max_count",
             "ff_log_factorial", "ff_log_coef", "ff_coef_power"))
}

# The data declarations can be shared with future NB/ZIP binomial families,
# but this function also installs Poisson-binomial Stan source, so is specific.
n_mixture_poisson_binomial_stanvars <- function(auxiliary) {
  source <- paste(readLines(system.file("stan", "n_mixture_poisson_binomial.stan",
    package = "flocker", mustWork = TRUE)), collapse = "\n")
  declarations <- c(
    ff_nunit = "int<lower=1> ff_nunit;",
    ff_ncoef = "int<lower=1> ff_ncoef;",
    ff_unit_offsets = "array[ff_nunit + 1] int<lower=1> ff_unit_offsets;",
    ff_coef_offsets = "array[ff_nunit + 1] int<lower=1> ff_coef_offsets;",
    ff_max_count = "array[ff_nunit] int<lower=0> ff_max_count;",
    ff_log_factorial = "vector[ff_nunit] ff_log_factorial;",
    ff_log_coef = "vector[ff_ncoef] ff_log_coef;",
    ff_coef_power = "vector[ff_ncoef] ff_coef_power;"
  )
  vars <- brms::stanvar(scode = source, block = "functions")
  for (name in names(declarations)) {
    value <- auxiliary[[name]]
    # rstan otherwise interprets length-one vectors as scalar data.
    if (!name %in% c("ff_nunit", "ff_ncoef")) value <- as.array(value)
    vars <- vars + brms::stanvar(value, name,
      scode = declarations[[name]], block = "data")
  }
  vars
}
