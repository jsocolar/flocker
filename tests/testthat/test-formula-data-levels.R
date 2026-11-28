test_that("transformed and combined formulas respect covariate levels", {
  fd <- make_flocker_data(
    matrix(c(1, 0, 0, 0, 1, 0), 3, 2),
    data.frame(g = factor(c("a", "a", "b")), x = 1:3),
    list(w = matrix(1:6, 3, 2)), type = "twolevel_single",
    level2_group = "g", quiet = TRUE
  )
  for (formula in list(~ log(w), ~ x:w, ~ (1 | w))) {
    expect_error(flocker_standata(f_occ = formula, f_det = ~ 1,
      f_meta = ~ 1, flocker_data = fd), "event covariates")
  }
  expect_error(flocker_standata(f_det = brms::bf(det ~ 1, occ ~ log(w), Omega ~ 1),
    flocker_data = fd), "event covariates")
  expect_error(flocker_standata(f_det = brms::bf(det ~ 1, occ ~ x, Omega ~ x),
    flocker_data = fd), "must be level-two covariates")
  expect_error(flocker_standata(f_det = brms::bf(det ~ 1, occ ~ x, Omega ~ g),
    flocker_data = fd), "cannot be used in f_meta")
  expect_no_error(flocker_standata(f_det = brms::bf(det ~ w, occ ~ x, Omega ~ 1),
    flocker_data = fd))

  for (i in 1:4) {
    formulas <- rep(list(NULL), 4)
    formulas[[i]] <- ~ log(w)
    expect_error(do.call(validate_unit_formula_variables,
      c(formulas, list(flocker_data = fd))), "event covariates")
  }
})
