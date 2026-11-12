test_that("single-season post-processing preserves names and singleton dimensions", {
  fit <- fixture_fit("single")
  skip_if(is.null(fit), "single-season fitted fixture unavailable")
  obs <- rbind(c(1, 0, 0), c(0, NA, NA))
  dimnames(obs) <- list(c("ash", "birch"), c("am", "noon", "pm"))
  unit <- data.frame(uc1 = c(0, 1), species = rep(fit$data$species[1], 2))
  ec <- matrix(seq_len(6) / 6, 2)
  for (n in c(2L, 1L)) {
    for (visits in c(3L, 1L)) {
      y <- obs[seq_len(n), seq_len(visits), drop = FALSE]
      fd <- make_flocker_data(y, unit[seq_len(n), , drop = FALSE],
                             list(ec1 = ec[seq_len(n), seq_len(visits), drop = FALSE]),
                             newdata_checks = TRUE, quiet = TRUE)
      expected_names <- c(dimnames(y), list("draw_1"))
      values <- fitted_flocker(fit, new_data = fd, draw_ids = 1)
      for (v in values) {
        expect_identical(dim(v), c(n, visits, 1L))
        expect_identical(dimnames(v), expected_names)
      }
      summary <- fitted_flocker(fit, new_data = fd, draw_ids = 1, summarise = TRUE)
      expect_identical(dim(summary$linpred_det), c(n, visits, 3L))
      expect_identical(dimnames(summary$linpred_det)[1:2], dimnames(y))
      uv <- fitted_flocker(fit, new_data = fd, components = "occ",
                           unit_level = TRUE, draw_ids = 1)$linpred_occ
      expect_identical(dimnames(uv), list(rownames(y), "draw_1"))
      for (condition in c(FALSE, TRUE)) {
        for (sample in c(FALSE, TRUE)) {
          z <- get_Z(fit, new_data = fd, draw_ids = 1,
                     history_condition = condition, sample = sample)
          expect_identical(dim(z), c(n, 1L))
          expect_identical(dimnames(z), list(rownames(y), "draw_1"))
          if (condition) expect_equal(unname(z[1, 1]), 1)
        }
        pred <- predict_flocker(fit, new_data = fd, draw_ids = 1,
                                history_condition = condition)
        expect_identical(dimnames(pred), expected_names)
        expect_identical(is.na(pred[, , 1, drop = FALSE]),
                         array(is.na(y), c(n, visits, 1L), expected_names))
      }
      ll <- log_lik_flocker(fit, new_data = fd, draw_ids = 1)
      expect_identical(dimnames(ll), list("draw_1", rownames(y)))
      expect_true(all(is.finite(ll)))
      mixed <- predict_flocker(fit, new_data = fd, draw_ids = 1, mixed = TRUE,
                               allow_new_levels = TRUE, sample_new_levels = "gaussian")
      expect_identical(dimnames(mixed), expected_names)
    }
  }
  before <- get_flocker_metadata(fit)
  training <- fitted_flocker(fit, draw_ids = 1)$linpred_det
  original_names <- before$original_dimnames
  for (i in seq_along(original_names)) {
    if (!is.null(original_names[[i]])) {
      expect_identical(dimnames(training)[[i]], original_names[[i]])
    }
  }
  expect_silent(predict_flocker(fit, draw_ids = 1, mixed = TRUE,
                               allow_new_levels = TRUE, sample_new_levels = "gaussian"))
  expect_identical(get_flocker_metadata(fit), before)
})

test_that("single-season schema checks cover all post-processing entry points", {
  fit <- fixture_fit("single")
  skip_if(is.null(fit), "single-season fitted fixture unavailable")
  fd <- structure(list(data = fit$data, metadata = get_flocker_metadata(fit)),
                  class = "flocker_data")
  functions <- list(fitted_flocker, get_Z, predict_flocker, log_lik_flocker)
  for (tag in list(NULL, 1L, 999L)) {
    bad <- fit
    attr(bad, "flocker_metadata")$schema_version <- tag
    bad_data <- fd
    bad_data$metadata$schema_version <- tag
    for (fun in functions) {
      expect_error(fun(bad), "Incompatible single-season packing schema")
      expect_error(fun(fit, new_data = bad_data), "Incompatible single-season packing schema")
    }
    expect_error(fitted_flocker(bad, new_data = fit$data), "Incompatible single-season packing schema")
    expect_error(get_Z(bad, new_data = fit$data, history_condition = FALSE),
                 "Incompatible single-season packing schema")
  }
  attr(bad, "flocker_metadata") <- NULL
  expect_error(fitted_flocker(bad, new_data = fit$data), "Incompatible single-season packing schema")
})
