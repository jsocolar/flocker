test_that("single-season source has a fixed input signature", {
  code <- make_occupancy_single_lpmf()
  expect_match(code, "array[] int unit_starts", fixed = TRUE)
  expect_no_match(code, "index_array", fixed = TRUE)
  expect_no_match(code, "vint4", fixed = TRUE)
})

test_that("augmented source delegates to the two-level likelihood", {
  code <- make_occupancy_augmented_lpmf()
  expect_match(code, make_occupancy_twolevel_single_lpmf(), fixed = TRUE)
  expect_match(code, "return occupancy_twolevel_single_lpmf(y |", fixed = TRUE)
})

test_that("two-level source reuses the single-season likelihood", {
  code <- make_occupancy_twolevel_single_lpmf()
  expect_match(code, make_occupancy_single_lpmf(), fixed = TRUE)
  expect_match(code, "occupancy_single_lpmf(y[group_start:group_end] |", fixed = TRUE)
})

test_that("no-argument functions have correct return types", {
  expect_is(make_occupancy_single_C_lpmf(), "character")
})

test_that("multiseason likelihood source has stable family entry points", {
  code <- make_occupancy_multi_lpmf()
  for (family in c("colex", "colex_eq", "autologistic", "autologistic_eq")) {
    expect_match(code, paste0("real occupancy_multi_", family, "_lpmf"), fixed = TRUE)
  }
})
