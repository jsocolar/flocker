test_that("fixed-input single-season source integrates with brms", {
  source <- paste(readLines(system.file(
    "stan", "occupancy_single.stan", package = "flocker", mustWork = TRUE
  )), collapse = "\n")
  family <- brms::custom_family(
    "occupancy_single", dpars = c("mu", "occ"),
    links = c("identity", "identity"), type = "int",
    vars = c("vint1", "vint2"), loop = FALSE
  )
  # Hand-packed units: three visits, one visit, and two visits. There is no
  # dependence on the existing formatter or its dimension-dependent indices.
  data <- data.frame(
    y = c(0L, 1L, 0L, 0L, 0L, 0L),
    unit_size = c(3L, -99L, -99L, 1L, 2L, -99L),
    Q = c(1L, 1L, 1L, 0L, 0L, 0L),
    x = c(10, 10, 10, 20, 30, 30),
    w = c(0.2, 0.5, 0.8, 0.4, 0.1, 0.7)
  )
  formula <- brms::bf(y | vint(unit_size, Q) ~ x + w, occ ~ x)
  stanvars <- brms::stanvar(scode = source, block = "functions")
  code <- brms::make_stancode(formula, data, family = family, stanvars = stanvars)
  expect_match(code, source, fixed = TRUE)
  expect_true(rstan::stanc(model_code = code)$status)
  # The same source also accepts a single unit with a single event.
  single <- data[4, , drop = FALSE]
  single$unit_size <- 1L
  single$Q <- as.integer(single$y > 0)
  code_single <- brms::make_stancode(
    brms::bf(y | vint(unit_size, Q) ~ 1, occ ~ 1), single,
    family = family, stanvars = stanvars
  )
  expect_match(code_single, source, fixed = TRUE)
  expect_true(rstan::stanc(model_code = code_single)$status)
  # When every unit has one event, starts fill the column with no padding.
  singletons <- data[1:3, , drop = FALSE]
  singletons$unit_size <- rep(1L, 3)
  singletons$Q <- as.integer(singletons$y > 0)
  code_singletons <- brms::make_stancode(
    brms::bf(y | vint(unit_size, Q) ~ 1, occ ~ 1), singletons,
    family = family, stanvars = stanvars
  )
  expect_true(rstan::stanc(model_code = code_singletons)$status)
})
