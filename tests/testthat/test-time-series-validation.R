test_that("time-series terms are rejected in every formula and model type", {
  for (type in flocker_data_output_types()) {
    for (name in c("f_occ", "f_det", "f_col", "f_ex", "f_auto", "f_meta")) {
      for (term in c("ar", "ma", "arma", "cosy", "unstr")) {
        args <- list(f_occ = ~ 1, f_det = ~ 1,
                     flocker_data = list(type = type), multiseason = NULL,
                     f_col = NULL, f_ex = NULL, multi_init = NULL,
                     f_auto = NULL, augmented = FALSE, threads = NULL,
                     f_meta = NULL)
        args[[name]] <- as.formula(paste0("~ ", term, "(time = t, gr = g)"))
        expect_error(do.call(validate_flock_params, args),
                     "time-series correlation terms")
      }
    }
  }
})

test_that("combined formulas and autocor specifications are inspected", {
  fd <- make_flocker_data(matrix(c(1, 0, 0, 1), 2), quiet = TRUE)
  formulas <- list(
    brms::bf(det ~ ar(time = t), occ ~ 1),
    brms::bf(det ~ 1, occ ~ brms::ar(time = t)),
    brms::bf(det ~ 1, occ ~ 1, autocor = ~ arma(time = t)),
    brms::bf(det ~ 1, occ ~ 1) + brms::acformula(~ ma(time = t)),
    brms::bf(det ~ 1, brms::lf(occ ~ 1, autocor = ~ ma(time = t)))
  )
  for (f in formulas) {
    expect_error(flocker_standata(f_det = f, flocker_data = fd),
                 "time-series correlation terms")
  }
  expect_error(flocker_standata(~ 1, ~ 1, fd, autocor = ~ ar(time = t)),
               "time-series correlation terms")
  expect_silent(validate_no_time_series_terms(~ ar + ma + arma + (1 | g)))
  expect_silent(validate_no_time_series_terms(~ car(M) + sar(M)))
  expect_silent(validate_no_time_series_terms(brms::bf(det ~ x, occ ~ x)))
  expect_silent(flocker_standata(~ 1, ~ 1, fd))
})
