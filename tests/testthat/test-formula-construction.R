test_that("formula responses are replaced without changing expressions or environments", {
  f <- local({
    passthrough <- function(x, ...) x
    ~ passthrough(`a b`, ~ 1) + (1 | group)
  })
  out <- set_formula_response(f, quote(ff_y | vint(ff_n_unit)))
  expect_identical(out[[2]], quote(ff_y | vint(ff_n_unit)))
  expect_identical(out[[3]], f[[2]])
  expect_identical(environment(out), environment(f))
  replaced <- set_formula_response(out, quote(det))
  expect_identical(replaced[[2]], quote(det))
  expect_identical(replaced[[3]], f[[2]])
  expect_identical(environment(replaced), environment(f))
})

test_that("all occupancy families retain local transformations and nested formulas", {
  # The nested ~ also exercises the former combined-formula string splitting.
  specs <- local({
    passthrough <- function(x, ...) x
    twice <- function(x) 2 * x
    list(
      det = ~ passthrough(x, ~ 1), process = ~ twice(x), meta = ~ twice(z),
      single = brms::bf(det ~ passthrough(x, ~ 1), occ ~ twice(x)),
      twolevel = brms::bf(det ~ passthrough(x, ~ 1), occ ~ twice(x), Omega ~ twice(z)),
      augmented = brms::bf(det ~ passthrough(x, ~ 1), occ ~ twice(x)),
      colex = brms::bf(det ~ passthrough(x, ~ 1), occ ~ twice(x),
                       colo ~ twice(x), ex ~ twice(x)),
      colex_eq = brms::bf(det ~ passthrough(x, ~ 1), colo ~ twice(x), ex ~ twice(x)),
      auto = brms::bf(det ~ passthrough(x, ~ 1), occ ~ twice(x),
                      colo ~ twice(x), autologistic ~ twice(x)),
      auto_eq = brms::bf(det ~ passthrough(x, ~ 1), colo ~ twice(x),
                         autologistic ~ twice(x))
    )
  })
  obs <- matrix(c(0, 1, 0, 1, 0, 1), 3)
  unit <- data.frame(x = 1:3, group = factor(c("a", "a", "b")))
  event <- list(w = matrix(1, 3, 2))
  multi <- array(rep(obs, 2), c(3, 2, 2))
  cases <- list(
    single = list(data = make_flocker_data(obs, unit, event, quiet = TRUE),
                  combined = "single", dpars = "occ"),
    single_C = list(data = make_flocker_data(obs, unit, quiet = TRUE),
                    combined = "single", dpars = "occ"),
    twolevel = list(data = make_flocker_data(obs, unit, event, type = "twolevel_single",
                                             level2_group = "group",
                                             level2_covs = data.frame(group = factor(c("a", "b")), z = 1:2),
                                             quiet = TRUE),
                    combined = "twolevel", dpars = c("occ", "Omega")),
    augmented = list(data = make_flocker_data(multi, unit, event, type = "augmented",
                                              n_aug = 2, quiet = TRUE),
                      combined = "augmented", dpars = "occ", augmented = TRUE)
  )
  fd <- make_flocker_data(multi, list(unit, unit),
                          list(w = array(1, dim(multi))), type = "multi", quiet = TRUE)
  for (kind in c("colex", "autologistic")) {
    for (init in c("explicit", "equilibrium")) {
      key <- paste(kind, init)
      cases[[key]] <- list(data = fd, multiseason = kind, multi_init = init,
        combined = paste0(if (kind == "colex") "colex" else "auto",
                           if (init == "equilibrium") "_eq" else ""),
        dpars = c(if (init == "explicit") "occ", "colo",
                  if (kind == "colex") "ex" else "autologistic"))
    }
  }
  arg_names <- c(occ = "f_occ", colo = "f_col", ex = "f_ex",
                 autologistic = "f_auto", Omega = "f_meta")
  for (case in cases) {
    args <- c(list(flocker_data = case$data, f_det = specs$det),
              case[intersect(names(case), c("multiseason", "multi_init", "augmented"))])
    for (dp in case$dpars) {
      args[[arg_names[[dp]]]] <- if (dp == "Omega") specs$meta else specs$process
    }
    actual <- do.call(flocker_standata, args)
    expect_equal(unname(actual$X[, 2]), as.numeric(case$data$data$x))
    for (dp in case$dpars) {
      expect_equal(unname(actual[[paste0("X_", dp)]][, 2]),
                    2 * as.numeric(case$data$data[[if (dp == "Omega") "z" else "x"]]))
    }
    args[names(args) %in% arg_names] <- NULL
    args$f_det <- specs[[case$combined]]
    expect_equal(do.call(flocker_standata, args), actual)
  }
})

test_that("combined formulas still require an unmodified det response", {
  d <- make_flocker_data(matrix(0, 2, 2), quiet = TRUE)
  expect_error(flocker_standata(f_det = brms::bf(y ~ 1, occ ~ 1), flocker_data = d),
               "beginning with `det ~`")
})

test_that("combined nonlinear specifications survive response replacement", {
  f <- local({
    passthrough <- function(x, ...) x
    brms::bf(det ~ passthrough(x, ~ 1), occ ~ a + b * x, a ~ 1, b ~ 1) +
      brms::set_nl(dpar = "occ")
  })
  d <- make_flocker_data(matrix(c(0, 1, 1, 0), 2),
                         data.frame(x = c(-1, 1)), quiet = TRUE)
  priors <- c(brms::set_prior("normal(0, 1)", nlpar = "a"),
              brms::set_prior("normal(0, 1)", nlpar = "b"))
  code <- flocker_stancode(f_det = f, flocker_data = d, prior = priors)
  expect_true(rstan::stanc(model_code = code)$status)
  expect_match(code, "nlp_a", fixed = TRUE)
  expect_match(code, "nlp_b", fixed = TRUE)
  f <- brms::bf(det ~ a + b * x, a ~ 1, b ~ 1, occ ~ 1, nl = TRUE) +
    brms::set_nl(FALSE, dpar = "occ")
  expect_identical(attributes(set_formula_response(f$formula, quote(ff_y))),
                    attributes(f$formula))
  code <- flocker_stancode(f_det = f, flocker_data = d, prior = priors)
  expect_true(rstan::stanc(model_code = code)$status)
})
