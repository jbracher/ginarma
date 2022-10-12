setwd("/home/johannes/Documents/underreporting/General_INARMA/R/Hermite")

# define true values for the three scenarios:
# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_psi <- c(0.5, 0.7, 0.9)
vals_phi <- c(0.5, 0.8, 0.9)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_xi <- 1 - vals_phi + vals_phi*vals_kappa

vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs


# summarize results:
# initialize lists
results_sim <- list()
removed_est_se <- list()
n_available <- matrix(ncol = 3, nrow = 3,
                      dimnames = list(vals_lgt, paste0("sc", 1:3)))

for(s in 1:3){ # loop over scenarios
  tau <- vals_tau[s]
  psi <- vals_psi[s]
  phi <- vals_phi[s]
  kappa <- vals_kappa[s]

  scenario_temp <- matrix(NA, nrow = 3, ncol = 17,
                          dimnames = list(c(250, 500, 1000),
                                          c("tau", "psi", "phi", "kappa", "T",
                                            "mean_tau", "sd_tau", "est_se_tau",
                                            "mean_psi", "sd_psi", "est_se_psi",
                                            "mean_phi", "sd_phi", "est_se_phi",
                                            "mean_kappa", "sd_kappa", "est_se_kappa")))
  scenario_temp[1, c("tau", "psi", "phi", "kappa")] <- c(tau, psi, phi, kappa)
  scenario_temp[, "T"] <- c(250, 500, 1000)
  removed_est_se_temp <- matrix(NA, ncol = 4, nrow = 3,
                                dimnames = list(c(250, 500, 1000), c("tau", "psi", "phi", "kappa")))

  # loop over lengths of time series:
  for(lgt in c(250, 500, 1000)){

    # get parameter estimates:
    pars <- read.csv(file = paste0("Results/he_estim_", tau, "_", psi, "_", phi, "_",
                                   kappa, "_", lgt, ".csv"))[, c("log_tau", "logit_psi", "logit_phi", "logit_kappa")]

    n_available_temp <- max(which(!is.na(pars[, 1])))
    pars <- pars[1:n_available_temp, ]
    n_available[as.character(lgt), s] <- n_available_temp

    # fill in means and standard deviations:
    scenario_temp[as.character(lgt), "mean_tau"] <- mean(exp(pars[, "log_tau"]), na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_psi"] <- mean(exp(pars[, "logit_psi"])/(1 + exp(pars[, "logit_psi"])), na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_phi"] <- mean(exp(pars[, "logit_phi"])/(1 + exp(pars[, "logit_phi"])), na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_kappa"] <- mean(exp(pars[, "logit_kappa"])/(1 + exp(pars[, "logit_kappa"])), na.rm = TRUE)

    scenario_temp[as.character(lgt), "sd_tau"] <- sd(exp(pars[, "log_tau"]), na.rm = TRUE)
    scenario_temp[as.character(lgt), "sd_psi"] <- sd(exp(pars[, "logit_psi"])/(1 + exp(pars[, "logit_psi"])), na.rm = TRUE)
    scenario_temp[as.character(lgt), "sd_phi"] <- sd(exp(pars[, "logit_phi"])/(1 + exp(pars[, "logit_phi"])), na.rm = TRUE)
    scenario_temp[as.character(lgt), "sd_kappa"] <- sd(exp(pars[, "logit_kappa"])/(1 + exp(pars[, "logit_kappa"])), na.rm = TRUE)

    # get estimated standard errors:
    ses <- read.csv(file = paste0("Results/he_ses_", tau, "_", psi, "_", phi, "_",
                                  kappa, "_", lgt, ".csv"))[1:n_available_temp, c("log_tau", "logit_psi", "logit_phi", "logit_kappa")]
    ses_delta <- NA*ses; colnames(ses_delta) <- c("tau", "psi", "phi", "kappa")
    # for tau:
    ses_delta[, "tau"] <- ses[, "log_tau"]*exp(pars[, "log_tau"])^2
    ses_delta[, "psi"] <- ses[, "logit_psi"]*exp(pars[, "logit_psi"])/(1 + exp(pars[, "logit_psi"]))^2
    # adapt for phi and kappa
    ses_delta[, "phi"] <- ses[, "logit_phi"]*exp(pars[, "logit_phi"])/(1 + exp(pars[, "logit_phi"]))^2
    ses_delta[, "kappa"] <- ses[, "logit_kappa"]*exp(pars[, "logit_kappa"])/(1 + exp(pars[, "logit_kappa"]))^2

    ses_cleaned <- ses_delta; ses_cleaned[ses_cleaned > 30] <- NA # set numerically instable ses to NA
    scenario_temp[as.character(lgt), c("est_se_tau", "est_se_psi", "est_se_phi", "est_se_kappa")] <-
      apply(ses_cleaned, 2, mean, na.rm = TRUE)

    # how many times was se estimation not possible or obviously instable?
    removed_est_se_temp[as.character(lgt), ] <- colSums(is.na(ses_cleaned))

  }
  # store results in list:
  results_sim[[s]] <- scenario_temp
  removed_est_se[[s]] <- removed_est_se_temp
}


# generate tex code:
for(s in 1:3){
  to_print <- matrix(as.character(format(results_sim[[s]],
                                          digits = 2,
                                          scientific = FALSE)),
                      nrow = 3)
  to_print[1, 1:4] <- c(vals_tau[s], vals_psi[s], vals_phi[s], vals_kappa[s])
  to_print[, 5] <- vals_lgt
  to_print[2:3, 1:4] <- ""

  # re-order columns:
  to_print <- to_print[, c(5, 1, 6:8,
                           2, 9:11,
                           3, 12:14,
                           4, 15:17)]

  library(xtable)
  write(
    print(xtable(to_print), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("../../Draft/table/sim_hermite_sc", s, ".tex")
  )
}


