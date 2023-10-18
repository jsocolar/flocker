## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----true-occupancy, results='hide', message=FALSE----------------------------
library(flocker); library(brms)
set.seed(3)

L <- 0.5
k <- .1
t0 <- -5
t <- seq(-15, 15, 1)
n_site_per_time <- 30
n_visit <- 3
det_prob <- .3

data <- data.frame(
  t = rep(t, n_site_per_time)
)

data$psi <- L/(1 + exp(-k*(t - t0)))
data$Z <- rbinom(nrow(data), 1, data$psi)
data$v1 <- data$Z * rbinom(nrow(data), 1, det_prob)
data$v2 <- data$Z * rbinom(nrow(data), 1, det_prob)
data$v3 <- data$Z * rbinom(nrow(data), 1, det_prob)

fd <- make_flocker_data(
  obs = as.matrix(data[,c("v1", "v2", "v3")]),
  unit_covs = data.frame(t = data[,c("t")]),
  event_covs <- list(dummy = matrix(rnorm(n_visit*nrow(data)), ncol = 3))
)


## ----fitting, results='hide', message=FALSE-----------------------------------
fit <- flock(f_det = brms::bf(
                 det ~ 1 + dummy,
                 occ ~ log(L/(1 + exp(-k*(t - t0)) - L)),
                 L ~ 1,
                 k ~ 1,
                 t0 ~ 1
               ) +
               brms::set_nl(dpar = "occ"),
             prior = 
               c(
                 prior(normal(0, 5), nlpar = "t0"),
                 prior(normal(0, 1), nlpar = "k"), 
                 prior(beta(1, 1), nlpar = "L", lb = 0, ub = 1)
                ),
             flocker_data = fd, 
             control = list(adapt_delta = 0.9),
             cores = 4)

## ----summary------------------------------------------------------------------
summary(fit)


## ----simulate-svc-------------------------------------------------------------
set.seed(1)
n <- 2000 # sample size
lscale <- 0.3 # square root of l of the gaussian kernel
sigma_gp <- 1 # sigma of the gaussian kernel
intercept <- 0 # occupancy logit-intercept
det_intercept <- -1 # detection logit-intercept
n_visit <- 4

# covariate data for the model
gp_data <- data.frame(
  x = rnorm(n), 
  y = rnorm(n),
  covariate = rnorm(n)
  )

# get distance matrix
dist.mat <- as.matrix(
  stats::dist(gp_data[,c("x", "y")])
  )

# get covariance matrix
cov.mat <- sigma_gp^2 * exp(- (dist.mat^2)/(2*lscale^2))

# simulate occupancy data
gp_data$coef <- mgcv::rmvn(1, rep(0, n), cov.mat)
gp_data$lp <- intercept + gp_data$coef * gp_data$covariate
gp_data$psi <- boot::inv.logit(gp_data$lp)
gp_data$Z <- rbinom(n, 1, gp_data$psi)

# simulate visit data
obs <- matrix(nrow = n, ncol = n_visit)
for(j in 1:n_visit){
  obs[,j] <- gp_data$Z * rbinom(n, 1, boot::inv.logit(det_intercept))
}


## ----fit-svc, results='hide', message=FALSE-----------------------------------
fd2 <- make_flocker_data(obs = obs, unit_covs = gp_data[, c("x", "y", "covariate")])
svc_mod <- flock(
  f_det = brms::bf(
                 det ~ 1,
                 occ ~ occint + g * covariate,
                 occint ~ 1,
                 g ~ 0 + gp(x, y, scale = FALSE, k = 20, c = 1.25)
               ) +
               brms::set_nl(dpar = "occ"),
  flocker_data = fd2,
  cores = 4
)


## ----summary-svc--------------------------------------------------------------
summary(svc_mod)


