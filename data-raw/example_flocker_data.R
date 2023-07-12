## code to prepare `example_flocker_data` dataset goes here
efd <- function(rep_constant = FALSE, seed = 123, n_pt = 30, n_sp = 30, 
                n_rep = 4) {
  if (!is.null(seed)) {set.seed(seed)}
  n_unit <- n_pt*n_sp
  
  backbone <- expand.grid(species = factor(paste0("sp_", 1:n_sp)), 
                          id_rep = 1:n_rep,
                          id_point = 1:n_pt)
  backbone$row_id <- seq_len(nrow(backbone))
  
  unit_covs <- within(unique(backbone[c("species", "id_point")]), {
    id_unit = 1:n_unit
    uc1 = stats::rnorm(n_pt)[id_point]
    uc2 = stats::rnorm(n_pt)[id_point]
    grp = sample(c(1:20), n_unit, replace = T)
    logit_psi = 0 + 1*uc1 - 1*uc2 + stats::rnorm(20)[grp]
    true_Z = stats::rbinom(n_unit, 1, boot::inv.logit(logit_psi))
  })
  
  bb <- merge(backbone, unit_covs)
  bb <- bb[order(bb$row_id), ]
  
  df_full <- within(bb, {
    id_point_rep <- interaction(id_point, id_rep)
    ec1 = stats::rnorm(n_pt * n_rep)[id_point_rep]
    ec2 = stats::rnorm(n_pt * n_rep)[id_point_rep]
    logit_theta = -.5 + .5*uc1 + stats::rnorm(20)[grp]
    if (!rep_constant) {
      logit_theta = logit_theta + .5*ec1 - 1*ec2
    }
    obs <- stats::rbinom(n_unit * n_rep, true_Z, boot::inv.logit(logit_theta))
  })
  
  # format data for return
  out <- list(obs = t(unstack(df_full[c("obs", "id_unit")], obs ~ id_unit)), 
              unit_covs = unit_covs[c("uc1", "uc2", "grp", "species")], 
              event_covs = list(
                ec1 = t(unstack(df_full[c("ec1", "id_unit")], ec1 ~ id_unit)), 
                ec2 = t(unstack(df_full[c("ec2", "id_unit")], ec2 ~ id_unit))
              ))
  rownames(out$obs) <- NULL
  rownames(out$event_covs$ec1) <- NULL
  rownames(out$event_covs$ec2) <- NULL
  return(out)
}

example_flocker_data <- efd()
example_flocker_data_constant <- efd(TRUE)

usethis::use_data(example_flocker_data)
usethis::use_data(example_flocker_data_constant)

