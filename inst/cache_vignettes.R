# This script creates .Rmd vignettes that include pre-computed R output so that
# whent they get built on CI it's a lightweight operation.

# This script should be run every time there is an update to the vignettes
# or an update to the package for which the vignettes are a desirable
# part of the test harness (but ideally the actual test harness should be 
# sufficient for testing, and we shouldn't rely on vignettes).

# Edits that touch only the text of the vignettes, and not the R computation,
# can be made with extreme care by changing both the .Rmd.orig files and the
# .Rmd files with computation cached.

old_wd <- getwd()

setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))

knitr::knit(
  "vignettes/augmented_models.Rmd.orig", 
  output = "vignettes/augmented_models.Rmd"
)

knitr::knit(
  "vignettes/flocker_tutorial.Rmd.orig", 
  output = "vignettes/flocker_tutorial.Rmd"
)

knitr::knit(
  "vignettes/nonlinear_models.Rmd.orig", 
  output = "vignettes/nonlinear_models.Rmd"
)

knitr::knit(
  "vignettes/sbc.Rmd.orig", 
  output = "vignettes/sbc.Rmd"
)

setwd(old_wd)
