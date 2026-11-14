test_that("fitted N-mixture components retain intrinsic dimensions and names", {
  obs <- rbind(oak = c(2, 0, NA), ash = c(0, 1, 3))
  colnames(obs) <- c("dawn", "noon", "dusk")
  d <- make_nmix_data(obs, quiet = TRUE)
  lp <- cbind(c(0, 0, log(3), log(3), log(3)),
              c(log(2), log(2), log(4), log(4), log(4)))
  format <- function(component, response = TRUE, summarise = FALSE) {
    format_nmix_fitted(lp, component, d$metadata, response, summarise,
                       c(0.05, 0.95), c(4, 2))
  }
  abund <- format("abund")
  expect_equal(unname(abund), rbind(c(1, 2), c(3, 4)))
  expect_identical(dimnames(abund), list(c("oak", "ash"), c("draw_4", "draw_2")))
  det <- format("det")
  expect_identical(dim(det), c(2L, 3L, 2L))
  expect_identical(dimnames(det), list(rownames(obs), colnames(obs), c("draw_4", "draw_2")))
  expect_true(all(is.na(det[1, 3, ])))
  expect_equal(unname(det[, 1, ]), stats::plogis(lp[c(1, 3), ]))
  expect_equal(unname(format("abund", response = FALSE)), lp[c(1, 3), ])
  summary <- format("abund", summarise = TRUE)
  expect_equal(unname(summary[, 1]), c(1.5, 3.5))
  expect_equal(unname(summary[1, ]), c(1.5, 1.05, 1.95))
  expect_identical(colnames(summary), c("mean", "Q5", "Q95"))
})

test_that("singleton dimensions and summary dimensions are retained", {
  d <- make_nmix_data(matrix(2, 1, 1), quiet = TRUE, newdata_checks = TRUE)
  for (component in c("abund", "det")) {
    out <- format_nmix_fitted(matrix(log(2), 1, 1), component, d$metadata,
                              TRUE, FALSE, c(0.05, 0.95), 7)
    expect_identical(dim(out), if (component == "abund") c(1L, 1L) else c(1L, 1L, 1L))
    summary <- format_nmix_fitted(matrix(log(2), 1, 1), component, d$metadata,
                                  TRUE, TRUE, c(0.05, 0.95), 7)
    expect_equal(as.vector(summary), rep(if (component == "abund") 2 else 2 / 3, 3))
    expect_identical(dim(summary), if (component == "abund") c(1L, 3L) else c(1L, 1L, 3L))
  }
  expect_error(fitted_nmix(d), "returned by nmix")
})
