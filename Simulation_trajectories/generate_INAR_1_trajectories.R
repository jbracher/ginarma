# Script generating the trajectories of Poisson INAR(1) models
library(inarma)
library(tidyverse)

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_kappa <- c(0.5, 0.6, 0.8)

vals_lgt <- c(250, 500, 1000)

n_sim <- 1e3  # Number of simulations

for (s in 1:3) {  # Loop over the scenarios

  # Grab the parameter values
  tau <- vals_tau[s]
  kappa <- vals_kappa[s]

  for (lgt in vals_lgt) {  # Loop over the lengths
    dat <- matrix(nrow = n_sim, ncol = lgt)

    for (k in 1:n_sim) {  # Loop over the iterations

      # Simulate an INARMA(1, 1) model
      set.seed(k)
      dat[k, ] <- sim_inarma(tau = tau, kappa = kappa, beta = 0,
                             family = "Poisson", offspring = "binomial",
                             lgt = lgt)$X

    }

    # Save in a file
    write_csv(
      as.data.frame(dat),
      file = paste0("Simulation_trajectories/sim_INAR1_s",
                    s, "_lgt_", lgt, "_pois.csv")
    )

  }
}