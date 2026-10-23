make_augmented_packing_fixture <- function(n_aug = 2L) {
  obs <- array(0, c(3, 3, 2), dimnames = list(c("alpha", "beta", "gamma"), c("dawn", "noon", "dusk"), c("wren", "thrush")))
  obs[1, 2:3, ] <- NA
  obs[3, 3, ] <- NA
  obs[2, 2, 1] <- 1
  obs[3, 1, 2] <- 1
  unit <- data.frame(uc1 = c(-1, 0, 1))
  event <- matrix(seq_len(9) / 10, 3)
  list(obs = obs, unit = unit, event = event,
       data = make_flocker_data(obs, unit, list(ec1 = event), type = "augmented",
                                n_aug = n_aug, quiet = TRUE))
}

test_that("augmented formatting is a two-level packing adapter", {
  x <- make_augmented_packing_fixture()
  fd <- x$data
  expect_named(fd, c("data", "metadata"))
  augmented_obs <- array(0, c(3, 3, 4))
  augmented_obs[, , 1:2] <- x$obs
  augmented_obs[, , 3:4] <- ifelse(is.na(x$obs), NA, 0)
  two <- make_flocker_data(expand_array_3D(augmented_obs),
    data.frame(group = factor(rep(1:4, each = 3)), uc1 = rep(x$unit$uc1, 4)),
    list(ec1 = stack_matrix(x$event, 4)), type = "twolevel_single",
    level2_group = "group", quiet = TRUE)
  columns <- setdiff(names(fd$data), "ff_species")
  expect_identical(fd$data[columns], two$data[columns])
  gp <- get_positions(fd)
  expect_equal(array(fd$data$ff_y[gp], dim(gp)), augmented_obs)
  expect_equal(dim(get_positions(fd, TRUE)), c(3, 4))
  expect_identical(fd$metadata$level2_group_names,
                   c("wren", "thrush", "pseudospecies1", "pseudospecies2"))
  starts <- fd$metadata$group_start
  expect_identical(fd$data$ff_group_Q[starts], c(1L, 1L, 0L, 0L))
  expect_identical(fd$data$ff_group_start[1:4], c(1L, 7L, 13L, 19L))
  expect_equal(fd$data$ff_n_unit[starts], rep(3L, 4))
  args <- list(f_occ = ~ uc1 + (1 | ff_species), f_det = ~ ec1,
               flocker_data = fd, augmented = TRUE)
  code <- do.call(flocker_stancode, args)
  expect_true(rstan::stanc(model_code = code)$status)
  sd <- do.call(flocker_standata, args)
  columns <- c("ff_n_group", "ff_group_start", "ff_n_unit", "ff_unit_start", "ff_Q", "ff_group_Q")
  for (i in 1:6) expect_equal(as.vector(sd[[paste0("vint", i)]]), fd$data[[columns[i]]])
  args$flocker_data <- make_augmented_packing_fixture(4)$data
  expect_identical(do.call(flocker_stancode, args), code)
  args$flocker_data$metadata$data_role <- "new_data"
  expect_error(do.call(flocker_standata, args), "new_data cannot")
  args$flocker_data$metadata$schema_version <- NULL
  expect_error(do.call(flocker_standata, args), "Incompatible augmented")
})

test_that("augmented post-processing preserves masks and both occupancy levels", {
  skip_on_cran()
  fd <- make_augmented_packing_fixture()$data
  fit <- example_flocker_model_aug
  for (condition in c(FALSE, TRUE)) {
    z <- get_Z(fit, new_data = fd, draw_ids = 1:2, history_condition = condition)
    expect_identical(dim(z$unit), c(3L, 4L, 2L))
    expect_identical(dimnames(z$unit)[1:2], fd$metadata$original_dimnames[c(1, 3)])
    expect_identical(dim(z$level2), c(4L, 2L))
    expect_identical(rownames(z$level2), fd$metadata$level2_group_names)
    p <- predict_flocker(fit, new_data = fd, draw_ids = 1:2, history_condition = condition)
    expect_identical(dim(p), c(3L, 3L, 4L, 2L))
    expect_identical(dimnames(p)[1:3], fd$metadata$original_dimnames)
    expect_equal(unname(is.na(p)), array(rep(is.na(get_positions(fd)), 2), dim(p)))
  }
  z <- get_Z(fit, new_data = fd, draw_ids = 1:2)
  expect_equal(unname(z$level2[1:2, ]), matrix(1, 2, 2))
  ll <- log_lik_flocker(fit, new_data = fd, draw_ids = 1:2)
  expect_identical(dim(ll), c(2L, 4L))
  expect_true(all(is.finite(ll)))
  expect_identical(colnames(ll), fd$metadata$level2_group_names)
})
