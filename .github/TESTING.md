# Running flocker tests

Load the checkout without loading test helpers or fixtures:

```r
pkgload::load_all(helpers = FALSE)
```

`tests/testthat/helper-fixtures.R` provides `fixture_data(name)` and
`fixture_fit(name)`. Data are simulated and formatted only on first request;
fits are compiled and sampled only on first request. A filtered data-only run
therefore does not construct any fitted fixtures:

```r
testthat::test_dir("tests/testthat", filter = "make_flocker_data")
```

Supported fixture names are `single`, `single_C`, `twolevel`, `aug`,
`multi_colex_ex`, `multi_colex_eq`, `multi_auto_ex`, and `multi_auto_eq`.
Request a fixture inside a test, after any skips. Getters also enforce the
existing platform/CRAN restrictions before looking up or constructing a fit.
Fixture generation preserves the caller's random-number state.

## Reusing fits locally

By default, fits are reused in memory for the current test run. For reuse
across sessions, explicitly configure a cache outside the package source:

```sh
FLOCKER_TEST_CACHE_DIR=/path/to/local/flocker-cache NOT_CRAN=true \
  Rscript -e 'pkgload::load_all(helpers=FALSE); testthat::test_dir("tests/testthat")'
```

Each family has its own RDS file. Cache signatures include the fixture version,
formatted input, formulas and sampling arguments, package R and Stan source,
OS/platform, R version, and relevant dependency versions. Missing, unreadable,
or incompatible entries are rebuilt only when that family is requested.
The old monolithic fixture cache is not reused. Cache files are development
artifacts and must not be committed or included in the source tarball.

For a cold integration run, unset `FLOCKER_TEST_CACHE_DIR` or point it at a new
empty directory. With `NOT_CRAN=true`, the full suite exercises all fitted
families. This lazy-fixture change does not change the existing CRAN test tier
or CI matrix; the smaller CRAN smoke-fit tier described in `V2_DESIGN.md` is
separate follow-up work.

SBC remains separate release validation. Ordinary tests and fixture generation
do not run or regenerate SBC vignettes.
