test_that("two-level groups contain standalone single-season datasets", {
  obs <- rbind(c(0, NA, NA), c(1, 0, 0), c(0, 0, NA), c(0, 0, 0))
  unit <- data.frame(g = factor(c("b", "a", "b", "a")), x = 1:4)
  event <- matrix(seq_along(obs), 4)
  fd <- make_flocker_data(obs, unit, list(w = event), type = "twolevel_single",
                          level2_group = "g", quiet = TRUE)
  expect_named(fd, c("data", "metadata"))
  expect_identical(fd$data$ff_group_start[1:2], c(1L, 7L))
  expect_identical(fd$data$ff_group_Q, c(rep(1L, 6), rep(0L, 3)))
  expect_equal(matrix(fd$data$ff_y[get_positions(fd)], 4), obs)
  for (g in 1:2) {
    units <- which(as.integer(unit$g) == g)
    visits <- seq_len(max(rowSums(!is.na(obs[units, , drop = FALSE]))))
    single <- make_flocker_data(obs[units, visits, drop = FALSE], unit[units, ],
      list(w = event[units, visits, drop = FALSE]), newdata_checks = TRUE, quiet = TRUE)
    rows <- if (g == 1) 1:6 else 7:9
    cols <- c("ff_y", "x", "w", "ff_n_unit", "ff_unit_start", "ff_Q")
    block <- fd$data[rows, cols]
    rownames(block) <- NULL
    expect_equal(block, single$data[cols])
  }
  args <- list(f_occ = ~ x, f_det = ~ w, f_meta = ~ 1, flocker_data = fd)
  code <- do.call(flocker_stancode, args)
  expect_true(rstan::stanc(model_code = code)$status)
  sd <- do.call(flocker_standata, args)
  cols <- c("ff_n_group", "ff_group_start", "ff_n_unit", "ff_unit_start", "ff_Q", "ff_group_Q")
  for (i in seq_along(cols)) expect_equal(as.vector(sd[[paste0("vint", i)]]), fd$data[[cols[i]]])
  larger <- rbind(obs, c(0, 0, 1))
  args$flocker_data <- make_flocker_data(larger,
    rbind(unit, data.frame(g = factor("b", levels = levels(unit$g)), x = 5)),
    list(w = matrix(1, 5, 3)), type = "twolevel_single", level2_group = "g", quiet = TRUE)
  expect_identical(do.call(flocker_stancode, args), code)
  args$flocker_data$metadata$data_role <- "new_data"
  expect_error(do.call(flocker_standata, args), "new_data cannot")
  args$flocker_data$metadata$schema_version <- NULL
  expect_error(do.call(flocker_standata, args), "Incompatible two-level")
  obs[1, ] <- NA
  expect_error(make_flocker_data(obs, unit, list(w = event), type = "twolevel_single",
                                level2_group = "g", quiet = TRUE), "first column")
})

test_that("two-level prediction handles a single unit and visit", {
  skip_on_cran()
  fd <- make_flocker_data(matrix(0, 1, 1), data.frame(group = factor("a")),
                          type = "twolevel_single", level2_group = "group",
                          newdata_checks = TRUE, quiet = TRUE)
  fit <- fixture_fit("twolevel")
  for (condition in c(FALSE, TRUE)) {
    z <- get_Z(fit, new_data = fd, draw_ids = 1:2, history_condition = condition)
    expect_identical(dim(z$unit), c(1L, 2L))
    expect_identical(dim(z$level2), c(1L, 2L))
    expect_true(all(is.finite(z$unit)))
    expect_true(all(z$unit <= z$level2))
    prediction <- predict_flocker(fit, new_data = fd, draw_ids = 1:2,
                                   history_condition = condition)
    expect_identical(dim(prediction), c(1L, 1L, 2L))
  }
  ll <- log_lik_flocker(fit, new_data = fd, draw_ids = 1:2)
  expect_identical(dim(ll), c(2L, 1L))
  expect_true(all(is.finite(ll)))
})
