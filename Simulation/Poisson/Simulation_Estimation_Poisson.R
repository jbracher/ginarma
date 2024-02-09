# choose simulation scenario and length of simulated time series
s <- 1
lgt <- 250

# scp johannes@130.60.71.234:/home/johannes/Documents/underreporting/General_INARMA/R/Poisson/Simulation_Estimation_Poisson.R Simulation_Estimation_Poisson.R
# scp johannes@130.60.71.234:/home/johannes/Documents/underreporting/General_INARMA/inarma_1.0.tar.gz inarma_1.0.tar.gz

continue <- FALSE # should existing results be read (to take up after error)

# setwd("/home/johannes/Documents/Projects/ingarch/Simulation/Poisson")
library(inarma)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs

# run simulation:
print(paste("Started scenario", s))

# grab true parameter values:
tau <- vals_tau[s]
beta <- vals_beta[s]
kappa <- vals_kappa[s]

# compute second-order properties (mentioned in manuscript):
sop <- sop_inarma(tau = tau, beta = beta, kappa = kappa, family = "Poisson")
sop

# initialize matrix to store results:
# get existing results if desired (if simulation somehow failed halfway)
if(continue){
  estim <- read.csv(file = paste0("Results/means_", tau, "_", beta, "_",
                                  kappa, "_", lgt, ".csv"))[, c("tau", "beta", "kappa", "mean_E1")]
  ses <- read.csv(file = paste0("Results/ses_", tau1, "_", phi, "_",
                                kappa, "_", lgt, ".csv"))[, c("tau", "beta", "kappa", "mean_E1")]
  estim_moments <- read.csv(file = paste0("Results/pois_estim_moments_", tau, "_", beta, "_",
                                          kappa, "_", lgt, ".csv"))[, c("tau", "beta", "kappa")]
  sim_data <- matrix(NA, ncol = lgt, nrow = n_sim)
}else{
  estim <- ses <- matrix(NA, nrow = n_sim, ncol = 4,
                         dimnames = list(NULL, c("tau", "beta", "kappa", "mean_E1")))
  estim_moments <- estim[, -ncol(estim)]
  sim_data <- matrix(NA, ncol = lgt, nrow = n_sim)
}

inds_to_run <- if(continue){
  (max(which(!is.na(estim[,"beta"]))) + 1):n_sim
}else{
  1:n_sim
}

# fit models:
for(i in inds_to_run){
  set.seed(i) # set seed
  # simulate all time series at once (wrote custom function for this)
  sim_temp <- sim_inarma(tau = tau,
                         beta = beta, kappa = kappa,
                         lgt = lgt, family = "Poisson")$X
  sim_data[i, ] <- sim_temp

  fit_temp <- NULL
  try({
    fit_temp <- fit_inarma(sim_temp, family = "Poisson", return_se = TRUE)
  })

  if(!is.null(fit_temp)){
    # store the parameter estimates
    estim[i, ] <- fit_temp$coefficients # store parameters
    # and store estimated standard errors:
    ses[i, ] <- fit_temp$se
  }

  # moment estimation:
  fit_temp_moments <- fit_inarma_moments(sim_temp, family = "Poisson")
  estim_moments[i, ] <- fit_temp_moments$coefficients # store parameters

  if(i%%10 == 0 | i < 5){
    # store:
    print(paste("Storing scenario", s, "after", i, "iterations..."))
    write.csv(sim_data, file = paste0("Results/pois_sim_", tau, "_", beta,
                                      "_", kappa, "_", lgt, ".csv"))
    write.csv(estim, file = paste0("Results/pois_estim_", tau, "_", beta,
                                   "_", kappa, "_", lgt, ".csv"))
    write.csv(ses, file = paste0("Results/pois_ses_", tau, "_", beta,
                                 "_", kappa, "_", lgt, ".csv"))
    write.csv(estim_moments, file = paste0("Results/pois_estim_moments_", tau, "_", beta,
                                           "_", kappa, "_", lgt, ".csv"))
  }
}
#   }
# }