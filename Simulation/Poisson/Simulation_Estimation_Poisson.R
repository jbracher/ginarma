# scp johannes@130.60.71.234:/home/johannes/Documents/underreporting/General_INARMA/R/Poisson/Simulation_Estimation_Poisson.R Simulation_Estimation_Poisson.R
# scp johannes@130.60.71.234:/home/johannes/Documents/underreporting/General_INARMA/inarma_1.0.tar.gz inarma_1.0.tar.gz

continue <- FALSE # should existing results be read (to take up after error)

# setwd("/home/johannes/Documents/underreporting/General_INARMA/R/Poisson")
library(inarma)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_phi <- c(0.5, 0.8, 0.9)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_xi <- 1 - vals_phi + vals_phi*vals_kappa

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs

vals_mu_I <- vals_sigma2_I <- vals_tau

vals_mu <- (vals_mu_I)/(1 - vals_kappa)
vals_sigma2 <- vals_sigma2_I/(1 - vals_kappa) -
  (vals_kappa*(1 + vals_xi - vals_phi*vals_kappa)*(vals_sigma2_I - vals_mu_I))/
  ((1 - vals_kappa)*(1 + vals_xi))
# identical to vals_mu for Poisson


# run simulation:
# for(lgt in vals_lgt){ # loop over lengths of time series
#   for(s in 1:3){ # loop over scenarios
print(paste("Started scenario", s))

# grab true parameter values:
tau <- vals_tau[s]
phi <- vals_phi[s]
kappa <- vals_kappa[s]


# initializ matrix to store results:
# get existing results if desired (if simulation somehow failed halfway)
if(continue){
  estim <- read.csv(file = paste0("Results/means_", tau1, "_", phi, "_",
                                  kappa, "_", lgt, ".csv"))[, c("log_tau", "logit_phi", "logit_kappa")]
  ses <- read.csv(file = paste0("Results/ses_", tau1, "_", phi, "_",
                                kappa, "_", lgt, ".csv"))[, c("log_tau", "logit_phi", "logit_kappa")]
  sim_data <- matrix(NA, ncol = lgt, nrow = n_sim)
}else{
  estim <- ses <- matrix(NA, nrow = n_sim, ncol = 3,
                         dimnames = list(NULL, c("log_tau", "logit_phi", "logit_kappa")))
  sim_data <- matrix(NA, ncol = lgt, nrow = n_sim)
}

inds_to_run <- if(continue){
  (max(which(!is.na(estim[,"phi"]))) + 1):n_sim
}else{
  1:n_sim
}

# fit models:
for(i in inds_to_run){
  set.seed(i) # set seed
  # simulate all time series at once (wrote custom function for this)
  sim_temp <- sim_inarma(tau = tau,
                         phi = phi, kappa = kappa,
                         lgt = lgt, family = "Poisson")$X
  sim_data[i, ] <- sim_temp

  fit_temp <- NULL
  try({
    fit_temp <- fit_inarma(sim_temp, family = "Poisson", return_se = TRUE, initialization = "stationary")
  })

  if(!is.null(fit_temp)){
    # store the parameter estimates
    estim[i, ] <- fit_temp$coefficients # store parameters
    # and store estimated standard errors:
    ses[i, ] <- sqrt(diag(solve(fit_temp$opt$hessian)))
  }

  if(i%%10 == 0 | i < 5){
    # store (commented out in order not to overwrite results):
    print(paste("Storing scenario", s, "after", i, "iterations..."))
    write.csv(sim_data, file = paste0("Results/Poisson_sim_", tau, "_", phi,
                                      "_", kappa, "_", lgt, ".csv"))
    write.csv(estim, file = paste0("Results/Poisson_estim_", tau, "_", phi,
                                   "_", kappa, "_", lgt, ".csv"))
    write.csv(ses, file = paste0("Results/Poisson_ses_", tau, "_", phi,
                                 "_", kappa, "_", lgt, ".csv"))
  }
}
#   }
# }