test_that("fixture caches are lazy, independent, and reused across sessions", {
  path <- tempfile()
  calls <- character()
  build <- function(name) { calls <<- c(calls, name); paste0("value-", name) }
  signature <- function(name) list(version = 1L, name = name)
  get <- new_fixture_cache(build, signature, path)
  expect_false(dir.exists(path))
  expect_length(calls, 0)
  expect_identical(get("a"), "value-a")
  expect_identical(get("a"), "value-a")
  expect_identical(calls, "a")
  expect_identical(list.files(path), "flocker-a.rds")
  second <- new_fixture_cache(build, signature, path)
  expect_identical(second("a"), "value-a")
  expect_identical(calls, "a")
  expect_identical(second("b"), "value-b")
  expect_identical(calls, c("a", "b"))
  changed <- new_fixture_cache(build, function(name) list(version = if (name == "a") 2L else 1L, name = name), path)
  expect_identical(changed("a"), "value-a")
  expect_identical(calls, c("a", "b", "a"))
  expect_identical(changed("b"), "value-b")
  expect_identical(calls, c("a", "b", "a"))
  writeLines("not an RDS", file.path(path, "flocker-b.rds"))
  repaired <- new_fixture_cache(build, signature, path)
  expect_identical(repaired("b"), "value-b")
  expect_identical(calls, c("a", "b", "a", "b"))
  unlink(path, recursive = TRUE)
})

test_that("memory-only caches do not require disk persistence", {
  calls <- 0L
  get <- new_fixture_cache(function(name) { calls <<- calls + 1L; name }, identity)
  expect_identical(get("single"), "single")
  expect_identical(get("single"), "single")
  expect_identical(calls, 1L)
})

test_that("all family data can be requested without fitting or changing RNG", {
  set.seed(321)
  seed <- .Random.seed
  for (name in fixture_names) {
    expect_silent(validate_packed_schema(fixture_data(name)))
    expect_identical(.Random.seed, seed)
  }
})

test_that("CRAN skips precede fixture construction", {
  withr::local_envvar(NOT_CRAN = "false")
  expect_condition(fixture_fit("twolevel"), class = "skip")
})
