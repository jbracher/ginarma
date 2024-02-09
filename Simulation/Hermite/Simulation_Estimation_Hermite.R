# choose simulation scenario and length of simulated time series
# this is usually done in a file calling this one, thus commented out.
# s <- 1
# lgt <- 250

# scp johannes@130.60.71.234:/home/johannes/Documents/underreporting/General_INARMA/R/Hermite/Simulation_Estimation_Hermite.R
# scp johannes@89.206.115.90:/home/johannes/Documents/underreporting/General_INARMA/inarma_1.0.tar.gz inarma_1.0.tar.gz

continue <- FALSE # should existing results be read (to take up after error)

# setwd("/home/johannes/Documents/underreporting/General_INARMA/R/Hermite")
library(inarma)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_psi <- c(0.5, 0.7, 0.9)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs




# run simulation:
print(paste("Started scenario", s))

# grab true parameter values:
tau <- vals_tau[s]
psi <- vals_psi[s]
beta <- vals_beta[s]
kappa <- vals_kappa[s]

# compute second-order properties (mentioned in manuscript):
sop <- sop_inarma(tau = tau, psi = psi, beta = beta, kappa = kappa, family = "NegBin")
sop


# initializ matrix to store results:
# get existing results if desired (if simulation somehow failed halfway)
if(continue){
  estim <- read.csv(file = paste0("Results/he_estim_", tau1, "_", tau2, "_", beta, "_",
                                  kappa, "_", lgt, ".csv"))[, c("tau", "beta", "kappa", "psi", "mean_E1")]
  estim_moments <- read.csv(file = paste0("Results/he_estim_moments_", tau, "_", psi, "_", beta, "_",
                                          kappa, "_", lgt, ".csv"))[, c("tau", "beta", "kappa", "psi")]
  ses <- read.csv(file = paste0("Results/he_ses_", tau1, "_", tau2, "_", beta, "_",
                                kappa, "_", lgt, ".csv"))[, c("tau", "beta", "kappa", "psi", "mean_E1")]
  sim_data <- matrix(NA, ncol = lgt, nrow = n_sim)
}else{
  estim <- ses <- matrix(NA, nrow = n_sim, ncol = 5,
                         dimnames = list(NULL, c("tau", "beta", "kappa", "psi", "mean_E1")))
  estim_moments <- matrix(NA, nrow = n_sim, ncol = 4,
                          dimnames = list(NULL, c("tau", "beta", "kappa", "psi")))
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
  sim_temp <- sim_inarma(tau = tau, psi = psi,
                         beta = beta, kappa = kappa,
                         lgt = lgt, family = "Hermite")$X
  sim_data[i, ] <- sim_temp
  
  fit_temp <- NULL
  try({
    fit_temp <- fit_inarma(sim_temp, family = "Hermite", return_se = TRUE)
    # very small values of overdispersion parameter indicate convergence issues.
    # Re-run with different starting values
    for(k in 1:4){
      if(fit_temp$coefficients_raw["logit_psi"] < -4){
        fit_temp_temp <- fit_inarma(sim_temp, family = "Hermite", return_se = TRUE,
                                    start = c("tau.Intercept" = c(0, 0.5, 1, 1.5)[k],
                                              "logit_phi" = c(0.3, 0.5, 0.7, 0.3)[k],
                                              "logit_kappa" = c(0.5, 0.6, 0.7, 0.4)[k],
                                              "logit_psi" = c(0, -1, 1, 0.5)[k]))
        if(fit_temp_temp$loglikelihood > fit_temp$loglikelihood){
          fit_temp <- fit_temp_temp
        }
      }
    }
  })

  if(!is.null(fit_temp)){
    # store the parameter estimates
    estim[i, ] <- fit_temp$coefficients # store parameters
    # and store estimated standard errors:
    ses[i, ] <- sqrt(diag(solve(fit_temp$opt$hessian)))
  }
  
  # moment estimation:
  fit_temp_moments <- fit_inarma_moments(sim_temp, family = "Hermite")
  estim_moments[i, ] <- fit_temp_moments$coefficients # store parameters
  
  if(i%%10 == 0 | i < 5){
    # store (commented out in order not to overwrite results):
    print(paste("Storing scenario", s, "after", i, "iterations..."))
    write.csv(sim_data, file = paste0("Results/he_sim_", tau, "_", psi, "_", beta,
                                      "_", kappa, "_", lgt, ".csv"))
    write.csv(estim, file = paste0("Results/he_estim_", tau, "_", psi, "_", beta,
                                  "_", kappa, "_", lgt, ".csv"))
    write.csv(ses, file = paste0("Results/he_ses_", tau, "_", psi, "_", beta,
                                 "_", kappa, "_", lgt, ".csv"))
    write.csv(estim_moments, file = paste0("Results/he_estim_moments_", tau, "_", psi, "_", beta,
                                           "_", kappa, "_", lgt, ".csv"))
  }
}