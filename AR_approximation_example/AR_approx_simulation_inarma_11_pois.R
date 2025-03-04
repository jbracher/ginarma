library(inarma)
library(tidyverse)

################################################################################
## Inference for the Poisson INARMA(1, 1) model using the approximation by an AR
## model
################################################################################

# Setup ------------------------------------------------------------------------

vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)

n_sim <- 1000 # number of simulation runs
lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000)  # lengths of simulated time series

# Run the loops ----------------------------------------------------------------

for (s in 1:3) {  # Loop over the scenarios

  # Grab the parameter values
  tau <- vals_tau[s]
  beta <- vals_beta[s]
  kappa <- vals_kappa[s]

  for (lgt in vals_lgt) {  # Loop over the series lengths

    # Allocate the result containers
    res_11 <- matrix(NA, nrow = n_sim, ncol = 6)
    colnames(res_11) <- c("tau", "kappa", "beta", "tau_se", "kappa_se", "beta_se")
    conv_11 <- rep(NA, n_sim)
    sims <- matrix(NA, nrow = n_sim, ncol = lgt)
    
    sim <- unname(as.matrix(read_csv(
      file = paste0("Simulation_trajectories/sim_INARMA11_s",
                    s, "_lgt_", lgt, "_pois.csv")
    )))

    for (k in 1:n_sim) {

      # Find the (approximate) maximum likelihood estimates
      fit <- fit_inarma_approx(sim$X, family = "Poisson", order = c(p = 1, q = 1),
                               lag_max = 10)

      # Save the results
      res_11[k, 1:3] <- fit$coefficients
      res_11[k, 4:6] <- fit$se
      conv_11[k] <- fit$convergence
      sims[k, ] <- sim$X
    }

    print(paste("Finished length:", lgt))

    results_tab <- as_tibble(res_11) %>% mutate(convergence = conv_11, length = lgt)
    sim_tab <- as_tibble(sims)

    # Write the results
    write_csv(results_tab, file = paste0("AR_approximation_example/Results/AR_INARMA11_s", s, "_lgt", lgt, "_pois.csv"))
  }
}
