test_that("single-level packing rejects empty units unless explicitly allowed", {
  obs <- rbind(c(1, 0), c(NA, NA), c(0, NA))
  event <- matrix(c(2, NA, 4, 6, NA, NA), 3)
  expect_error(pack_unit_block(obs, event_covs = list(w = event)), "Empty units")
  expect_error(pack_nested_blocks(obs, NULL, list(w = event), c(1, 1, 2)), "Empty units")
  expect_error(pack_unit_block(obs, event_covs = list(w = event), allow_empty = TRUE),
               "replacement value")
  packed <- pack_unit_block(obs, event_covs = list(w = event), allow_empty = TRUE,
                             placeholder_values = list(w = 12))
  expect_identical(packed$data$ff_y, c(1L, 0L, -99L, 0L))
  expect_equal(packed$data$w, c(2, 6, 12, 4))
  expect_equal(matrix(packed$data$ff_y[packed$original_to_packed], 3), obs)
})

test_that("multiseason placeholder imputation stays global across series", {
  obs <- array(0, c(2, 2, 3))
  obs[, , 2] <- NA
  event <- array(NA_real_, dim(obs))
  event[1, , c(1, 3)] <- 2
  event[2, , c(1, 3)] <- 10
  observer <- array(NA_character_, dim(obs))
  observer[1, , c(1, 3)] <- "z"
  observer[2, , c(1, 3)] <- "a"
  fd <- make_flocker_data(obs, event_covs = list(w = event, observer = observer),
                          type = "multi", quiet = TRUE, newdata_checks = TRUE)
  placeholders <- fd$metadata$placeholder
  expect_equal(fd$data$w[placeholders], c(6, 6))
  expect_equal(as.character(fd$data$observer[placeholders]), c("a", "a"))
  expect_identical(levels(fd$data$observer), c("a", "z"))

  # Each extracted series is exactly a one-level block with shared replacements.
  for (s in 1:2) {
    one <- pack_unit_block(t(obs[s, , ]),
      event_covs = list(w = t(event[s, , ]), observer = t(observer[s, , ])),
      allow_empty = TRUE, placeholder_values = list(w = 6, observer = "a"))
    start <- fd$metadata$series_start[s]
    rows <- seq.int(start, length.out = nrow(one$data))
    block <- fd$data[rows, names(one$data)]
    block$observer <- as.character(block$observer)
    rownames(block) <- NULL
    expect_equal(block, one$data)
  }
})
