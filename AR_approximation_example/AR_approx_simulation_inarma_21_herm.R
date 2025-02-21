library(inarma)
library(tidyverse)

##############################################################################
## Inference for the Hermite INARMA(2, 1) model using the approximation by an
## AR model
##############################################################################

# Setup ------------------------------------------------------------------------

lag_max <- 12  # The maximum lag of the approximating AR model
vals_lgt <- c(250, 500, 1000, 2000)

# Run the loops ----------------------------------------------------------------

n_sim <- 1000  # Number of iterations

for (s in 1:2) {
  for (lgt in vals_lgt) {  # Loop over the series lengths
    
    # Allocate the result containers
    res_21 <- matrix(NA, nrow = n_sim, ncol = 10)
    colnames(res_21) <- c("tau", "kappa1", "kappa2", "beta", "psi", "tau_se", "kappa1_se",
                          "kappa2_se", "beta_se", "psi_se")
    conv_21 <- rep(NA, n_sim)
    
    # Load the simulation trajectories
    dat <- unname(as.matrix(read.csv(
      file = paste0("Simulation_trajectories/sim_INARMA21_s",
                    s, "_lgt_", lgt, "_herm.csv"))))
    
    for (k in 1:n_sim) {
      
      # Find the (approximate) maximum likelihood estimates
      fit <- fit_inarma_approx(dat[k, ], family = "Hermite", order = c(p = 2, q = 1),
                               lag_max = lag_max, control_optim = list(maxit = 1200))
      
      # Save the results
      res_21[k, 1:5] <- fit$coefficients
      res_21[k, 6:10] <- fit$se
      conv_21[k] <- fit$convergence
      if (k %% 100 == 0) {
        print(paste0("Finished iteration ", k))
      }
    }
    
    print(paste("Finished length:", lgt))
    
    results_tab <- as_tibble(res_21) %>% mutate(convergence = conv_21, length = lgt)
    
    # Write the results
    write_csv(results_tab, file = paste0("AR_approximation_example/Results/AR_INARMA21_s", s, "_lgt_", lgt, "_herm.csv"))
  }
}