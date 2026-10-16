library(testthat)
library(flocker)

windows_cran <- identical(.Platform$OS.type, "windows") &&
  !identical(Sys.getenv("NOT_CRAN"), "true")

if (windows_cran) {
  message(
    "Skipping tests that require runtime Stan compilation on CRAN Windows"
  )
  test_check(
    "flocker",
    filter = paste0(
      "(exact-likelihoods|flock|make_flocker_data|",
      "make_occupancy_lpmf_functions|occupancy_families|",
      "simulate_flocker_data|utils)$"
    )
  )
} else {
  test_check("flocker")
}
