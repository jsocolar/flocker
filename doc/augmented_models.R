## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data-augmented, echo = FALSE---------------------------------------------
library(flocker)
d <- simulate_flocker_data(
    augmented = TRUE
    )
fd = make_flocker_data(
  d$obs, d$unit_covs, d$event_covs,
  type = "augmented", n_aug = 100
  )
fm <- flock(
  f_occ = ~ (1 | ff_species),
  f_det = ~ uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
  flocker_data = fd,
  augmented = TRUE,
  cores = 4
  )


## ----summary------------------------------------------------------------------
summary(fm)


