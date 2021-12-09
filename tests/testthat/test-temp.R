test_that("temporary test", {
  library(brms)
  df <- data.frame(x = rnorm(100), 
                   y = rnorm(100))
  f <- brm(y ~ x, data = df)
})