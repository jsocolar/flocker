test_that("multiseason blocks nest series, units, and visits", {
  x <- make_ragged_multi_fixture(global_trailing_season = TRUE)
  dimnames(x$obs) <- list(letters[1:3], paste0("v", 1:4), paste0("t", 1:6))
  # Values at absent visits must not influence numeric imputation.
  x$event[is.na(x$obs)] <- 10000
  x$event[1, 1, 2] <- NA
  fd <- suppressWarnings(make_flocker_data(x$obs, x$unit_covs,
                                           list(w = x$event), type = "multi", quiet = TRUE))
  m <- fd$metadata
  expect_named(fd, c("data", "metadata"))
  expect_identical(m$schema_version, 2L)
  expect_identical(m$data_role, "model_data")
  expect_identical(names(fd$data), c("ff_y", "uc1", "w", "ff_n_series",
                                    "ff_series_start", "ff_n_unit", "ff_unit_start", "ff_Q"))
  expect_identical(fd$data$ff_n_series, c(3L, rep(-99L, nrow(fd$data) - 1)))
  expect_identical(fd$data$ff_series_start[1:3], c(1L, 12L, 21L))
  expect_identical(m$series_start, c(1L, 12L, 21L))
  expect_identical(fd$data$ff_n_unit[c(1, 12, 21)], c(5L, 4L, 5L))
  expect_true(all(fd$data$ff_n_unit[-c(1, 12, 21)] == -99L))
  starts <- list(c(1L, 5L, 6L, 8L, 9L), c(1L, 4L, 5L, 6L),
                 c(1L, 3L, 7L, 10L, 11L))
  for (s in 1:3) {
    rows <- seq.int(m$series_start[s], c(11L, 20L, 31L)[s])
    block <- fd$data[rows, ]
    units <- block$ff_n_unit[1]
    expect_identical(block$ff_unit_start[seq_len(units)], starts[[s]])
    expect_true(all(block$ff_unit_start[-seq_len(units)] == -99L))
    for (u in seq_len(units)) {
      unit_end <- if (u < units) starts[[s]][u + 1L] - 1L else nrow(block)
      y <- block$ff_y[seq.int(starts[[s]][u], unit_end)]
      expect_identical(block$ff_Q[seq.int(starts[[s]][u], unit_end)],
                       rep(as.integer(any(y == 1)), length(y)))
    }
  }
  expect_equal(array(fd$data$ff_y[get_positions(fd)], dim(x$obs), dimnames(x$obs)), x$obs)
  expect_equal(unname(matrix(fd$data$uc1[get_positions(fd, TRUE)], 3)), x$expected_unit)
  expect_identical(dimnames(get_positions(fd)), dimnames(x$obs))
  expect_identical(dimnames(get_positions(fd, TRUE)), dimnames(x$obs)[c(1, 3)])
  expect_equal(fd$data$w[get_positions(fd, TRUE)[1, 2]], mean(x$event[!is.na(x$obs)]))
  expect_true(all(fd$data$w[m$placeholder][-1] == 10000))
  expect_true(all(fd$data$ff_y[m$placeholder] == -99L))
  expect_equal(sum(!m$placeholder), sum(!is.na(x$obs)))
})

test_that("factor placeholders retain reference levels and supplied values", {
  x <- make_ragged_multi_fixture()
  values <- factor(rep(c("second", "reference"), length.out = length(x$obs)),
                   levels = c("reference", "second"))
  dim(values) <- dim(x$obs)
  values[is.na(x$obs)] <- NA
  values[2, 1, 3] <- "second"
  fd <- make_flocker_data(x$obs, x$unit_covs, list(observer = values), type = "multi", quiet = TRUE)
  expect_true(is.factor(fd$data$observer))
  expect_identical(levels(fd$data$observer), levels(values))
  expect_identical(as.character(fd$data$observer[get_positions(fd, TRUE)[1, 2]]), "reference")
  expect_identical(as.character(fd$data$observer[get_positions(fd, TRUE)[2, 3]]), "second")
  values[1, 1, 1] <- NA
  expect_error(make_flocker_data(x$obs, x$unit_covs, list(observer = values), type = "multi", quiet = TRUE),
               "where the response is not missing")
})

test_that("all four multiseason families use stable programs and aligned data", {
  for (model in c("colex", "autologistic")) for (init in c("explicit", "equilibrium")) {
    args <- list(f_det = ~ ec1, f_col = ~ uc1, multiseason = model, multi_init = init)
    if (init == "explicit") args$f_occ <- ~ uc1
    if (model == "colex") args$f_ex <- ~ uc1 else args$f_auto <- ~ uc1
    args$flocker_data <- format_ragged_multi_fixture()
    code <- do.call(flocker_stancode, args)
    expect_true(rstan::stanc(model_code = code)$status)
    sd <- do.call(flocker_standata, args)
    d <- args$flocker_data$data
    for (i in 1:5) expect_equal(as.vector(sd[[paste0("vint", i)]]), d[[i + 3]])
    expect_equal(as.vector(sd$Y), d$ff_y)
    expect_equal(unname(sd$X_colo[, "uc1"]), d$uc1)
    args$flocker_data <- make_flocker_data(array(0, c(1, 2, 2)),
      list(data.frame(uc1 = 0), data.frame(uc1 = 1)), list(ec1 = array(1:4, c(1, 2, 2))),
      type = "multi", quiet = TRUE)
    expect_identical(do.call(flocker_stancode, args), code)
    old <- args$flocker_data
    args$flocker_data$metadata$schema_version <- NULL
    expect_error(do.call(flocker_standata, args), "Incompatible multiseason packing schema")
    expect_error(get_positions(args$flocker_data), "Incompatible multiseason packing schema")
    args$flocker_data <- old
    args$flocker_data$metadata$data_role <- "new_data"
    expect_error(do.call(flocker_standata, args), "new_data cannot")
  }
})

test_that("nested packing also supports unordered group membership", {
  obs <- rbind(c(0, 1), c(1, NA), c(0, 0))
  packed <- pack_nested_blocks(obs, data.frame(x = 1:3), list(w = obs + 1), c(2L, 1L, 2L))
  expect_identical(packed$unit_order, c(2L, 1L, 3L))
  expect_identical(packed$group_start, c(1L, 2L))
  expect_equal(matrix(packed$data$ff_y[packed$original_to_packed], 3), obs)
})

test_that("multiseason linear indexing handles three visits and singleton dimensions", {
  for (n_series in c(1L, 3L)) for (n_visit in 1:3) for (n_season in 1:3) {
    dims <- c(n_series, n_visit, n_season)
    obs <- array(rep(0:1, length.out = prod(dims)), dims)
    event <- array(seq_len(prod(dims)), dims)
    fd <- make_flocker_data(obs, lapply(seq_len(n_season), function(t) {
      data.frame(x = seq_len(n_series) + 10 * t)
    }), list(w = event), type = "multi", newdata_checks = TRUE, quiet = TRUE)
    expect_equal(array(fd$data$ff_y[get_positions(fd)], dims), obs)
    expect_equal(array(fd$data$w[get_positions(fd)], dims), event)
    expect_equal(matrix(fd$data$x[get_positions(fd, TRUE)], n_series),
                 outer(seq_len(n_series), 10 * seq_len(n_season), `+`))
  }
})

test_that("multiseason post-processing preserves singleton axes and input names", {
  skip_on_cran()
  fits <- list(example_flocker_model_multi_colex_ex, example_flocker_model_multi_colex_eq,
               example_flocker_model_multi_auto_ex, example_flocker_model_multi_auto_eq)
  for (n_season in 1:2) {
    obs <- array(0, c(1, 1, n_season),
                 dimnames = list("site", "visit", paste0("season", seq_len(n_season))))
    fd <- make_flocker_data(obs, rep(list(data.frame(uc1 = 0)), n_season),
                            list(ec1 = obs), type = "multi", newdata_checks = TRUE, quiet = TRUE)
    for (fit in fits) {
      for (condition in c(FALSE, TRUE)) {
        z <- get_Z(fit, draw_ids = 1:2, new_data = fd, history_condition = condition)
        expect_identical(dim(z), c(1L, as.integer(n_season), 2L))
        expect_identical(dimnames(z)[1:2], dimnames(obs)[c(1, 3)])
        expect_true(all(is.finite(z)))
        p <- predict_flocker(fit, draw_ids = 1:2, new_data = fd, history_condition = condition)
        expect_identical(dim(p), c(dim(obs), 2L))
        expect_identical(dimnames(p)[1:3], dimnames(obs))
      }
      ll <- log_lik_flocker(fit, draw_ids = 1:2, new_data = fd)
      expect_identical(dim(ll), c(2L, 1L))
      expect_identical(colnames(ll), "site")
      expect_true(all(is.finite(ll)))
    }
  }
})

test_that("a one-season forward-backward calculation reduces to Bayes' rule", {
  expect_equal(forward_backward_algorithm(1, 0.6, 0.3, 0.2, 0.4),
               0.3 * 0.6 / (0.7 + 0.3 * 0.6))
  expect_equal(forward_backward_algorithm(0, 0.4, 0.3, 0.2, 0.4), 1)
})

test_that("a series unit directory is unchanged when packed on its own", {
  x <- make_ragged_multi_fixture()
  together <- format_ragged_multi_fixture()
  starts <- together$data$ff_series_start[1:3]
  ends <- c(starts[-1] - 1L, nrow(together$data))
  for (s in 1:3) {
    alone <- suppressWarnings(make_flocker_data(
      x$obs[s, , , drop = FALSE],
      lapply(x$unit_covs, function(d) d[s, , drop = FALSE]),
      list(ec1 = x$event[s, , , drop = FALSE]), type = "multi", quiet = TRUE))
    columns <- c("ff_y", "uc1", "ff_n_unit", "ff_unit_start", "ff_Q")
    block <- together$data[seq.int(starts[s], ends[s]), columns]
    rownames(block) <- NULL
    expect_identical(block, alone$data[columns])
  }
})
