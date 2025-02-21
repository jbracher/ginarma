library(inarma)
library(tidyverse)

##############################################################################
## Inference for the NegBin INARMA(1, 1) model using the approximation by an
## AR model
##############################################################################

# Setup ------------------------------------------------------------------------

n_sim <- 1000 # number of simulation runs
lag_max <- 10  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000)  # lengths of simulated time series

# Run the loops ----------------------------------------------------------------

for (s in 1:3) {  # Loop over the scenarios

  for (lgt in vals_lgt) {  # Loop over the series lengths

    # Allocate the result containers
    res_11 <- matrix(NA, nrow = n_sim, ncol = 8)
    colnames(res_11) <- c("tau", "kappa", "beta", "psi", "tau_se", "kappa_se",
                          "beta_se", "psi_se")
    conv_11 <- rep(NA, n_sim)
    sims <- matrix(NA, nrow = n_sim, ncol = lgt)
    
    sim <- unname(as.matrix(read_csv(
      file = paste0("Simulation_trajectories/sim_INARMA11_s",
                    s, "_lgt_", lgt, "_nbin.csv")
    )))

    for (k in 1:n_sim) {

      # Find the (approximate) maximum likelihood estimates
      fit <- fit_inarma_approx(sim$X, family = "NegBin", order = c(p = 1, q = 1),
                               lag_max = 10)

      # Save the results
      res_11[k, 1:4] <- fit$coefficients
      res_11[k, 5:8] <- fit$se
      conv_11[k] <- fit$convergence
      sims[k, ] <- sim$X
    }

    print(paste("Finished length:", lgt))

    results_tab <- as_tibble(res_11) %>% mutate(convergence = conv_11, length = lgt)
    sim_tab <- as_tibble(sims)

    # Write the results
    write_csv(results_tab, file = paste0("AR_approximation_example/Results/AR_INARMA11_s", s, "_lgt", lgt, "_nbin.csv"))
  }
}
