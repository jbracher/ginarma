setwd("/home/johannes/Documents/Projects/ingarch/Simulation/NegBin")

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

  scenario_temp <- matrix(NA, nrow = 3, ncol = 25,
                          dimnames = list(c(250, 500, 1000),
                                          c("tau", "psi", "phi", "kappa", "T",
                                            "mean_tau", "se_tau", "est_se_tau", "mean_moment_tau", "se_moment_tau",
                                            "mean_psi", "se_psi", "est_se_psi", "mean_moment_psi", "se_moment_psi",
                                            "mean_phi", "se_phi", "est_se_phi", "mean_moment_phi", "se_moment_phi",
                                            "mean_kappa", "se_kappa", "est_se_kappa", "mean_moment_kappa", "se_moment_kappa")))
  scenario_temp[1, c("tau", "psi", "phi", "kappa")] <- c(tau, psi, phi, kappa)
  scenario_temp[, "T"] <- c(250, 500, 1000)
  removed_est_se_temp <- matrix(NA, ncol = 4, nrow = 3,
                                dimnames = list(c(250, 500, 1000), c("tau", "psi", "phi", "kappa")))

  # loop over lengths of time series:
  for(lgt in c(250, 500, 1000)){

    # get parameter estimates:
    pars <- read.csv(file = paste0("Results/nb_estim_", tau, "_", psi, "_", phi, "_",
                                   kappa, "_", lgt, ".csv"))[, c("log_tau", "log_psi", "logit_phi", "logit_kappa")]

    n_available_temp <- max(which(!is.na(pars[, 1])))
    pars <- pars[1:n_available_temp, ]
    n_available[as.character(lgt), s] <- n_available_temp

    # fill in means and standard deviations:
    scenario_temp[as.character(lgt), "mean_tau"] <- mean(exp(pars[, "log_tau"]), na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_psi"] <- mean(exp(pars[, "log_psi"]), na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_phi"] <- mean(exp(pars[, "logit_phi"])/(1 + exp(pars[, "logit_phi"])), na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_kappa"] <- mean(exp(pars[, "logit_kappa"])/(1 + exp(pars[, "logit_kappa"])), na.rm = TRUE)

    scenario_temp[as.character(lgt), "se_tau"] <- sd(exp(pars[, "log_tau"]), na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_psi"] <- sd(exp(pars[, "log_psi"]), na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_phi"] <- sd(exp(pars[, "logit_phi"])/(1 + exp(pars[, "logit_phi"])), na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_kappa"] <- sd(exp(pars[, "logit_kappa"])/(1 + exp(pars[, "logit_kappa"])), na.rm = TRUE)

    # get estimated standard errors:
    ses <- read.csv(file = paste0("Results/nb_ses_", tau, "_", psi, "_", phi, "_",
                                  kappa, "_", lgt, ".csv"))[1:n_available_temp, c("log_tau", "log_psi", "logit_phi", "logit_kappa")]
    ses_delta <- NA*ses; colnames(ses_delta) <- c("tau", "psi", "phi", "kappa")
    # for tau:
    ses_delta[, "tau"] <- ses[, "log_tau"]*exp(pars[, "log_tau"])^2
    ses_delta[, "psi"] <- ses[, "log_psi"]*exp(pars[, "log_psi"])^2
    # adapt for phi and kappa
    ses_delta[, "phi"] <- ses[, "logit_phi"]*exp(pars[, "logit_phi"])/(1 + exp(pars[, "logit_phi"]))^2
    ses_delta[, "kappa"] <- ses[, "logit_kappa"]*exp(pars[, "logit_kappa"])/(1 + exp(pars[, "logit_kappa"]))^2

    ses_cleaned <- ses_delta; ses_cleaned[ses_cleaned > 30] <- NA # set numerically instable ses to NA
    scenario_temp[as.character(lgt), c("est_se_tau", "est_se_psi", "est_se_phi", "est_se_kappa")] <-
      apply(ses_cleaned, 2, mean, na.rm = TRUE)

    # how many times was se estimation not possible or obviously instable?
    removed_est_se_temp[as.character(lgt), ] <- colSums(is.na(ses_cleaned))

    # get moment estimates:
    pars_moment <- read.csv(file = paste0("Results/nb_estim_moments_", tau, "_", psi, "_", phi, "_",
                                   kappa, "_", lgt, ".csv"))[, c("tau", "psi", "phi", "kappa")]

    scenario_temp[as.character(lgt), "mean_moment_tau"] <- mean(pars_moment[, "tau"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_moment_psi"] <- mean(pars_moment[, "psi"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_moment_phi"] <- mean(pars_moment[, "phi"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "mean_moment_kappa"] <- mean(pars_moment[, "kappa"], na.rm = TRUE)

    scenario_temp[as.character(lgt), "se_moment_tau"] <- sd(pars_moment[, "tau"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_moment_psi"] <- sd(pars_moment[, "psi"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_moment_phi"] <- sd(pars_moment[, "phi"], na.rm = TRUE)
    scenario_temp[as.character(lgt), "se_moment_kappa"] <- sd(pars_moment[, "kappa"], na.rm = TRUE)

  }
  # the package uses a different parameterization to the text with phi = 1 - beta
  # adapt the outputs accordingly:
  colnames(scenario_temp)[colnames(scenario_temp) == "phi"] <- "beta"
  scenario_temp[, "beta"] <- 1 - scenario_temp[, "beta"]

  scenario_temp[, "mean_phi"] <- 1 - scenario_temp[, "mean_phi"]
  colnames(scenario_temp)[colnames(scenario_temp) == "mean_phi"] <- "mean_beta"
  # for the standard errors only need to change labelling (as beta = 1 - phi)
  colnames(scenario_temp)[colnames(scenario_temp) == "se_phi"] <- "se_beta"
  colnames(scenario_temp)[colnames(scenario_temp) == "est_se_phi"] <- "est_se_beta"

  # fill in moment estimators:
  scenario_temp[, "mean_moment_phi"] <- 1 - scenario_temp[, "mean_moment_phi"]
  colnames(scenario_temp)[colnames(scenario_temp) == "mean_moment_phi"] <- "mean_moment_beta"
  colnames(scenario_temp)[colnames(scenario_temp) == "se_moment_phi"] <- "se_moment_beta"

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
  colnames(to_print) <- colnames(results_sim[[s]])
  # to_print[1, 1:4] <- c(vals_tau[s], vals_psi[s], 1 - vals_phi[s], vals_kappa[s]) # beta = 1 - phi
  to_print[, 5] <- vals_lgt
  to_print[2:3, 1:4] <- ""

  # re-order columns for table on ML estimate:
  to_print_ml <- to_print[, c(5, 1, 6:8,
                           2, 11:13,
                           3, 16:18,
                           4, 21:23)]
  to_print_ml <- to_print[, c("T",
                              "tau", "mean_tau", "se_tau", "est_se_tau",
                              "psi", "mean_psi", "se_psi", "est_se_psi",
                              "beta", "mean_beta", "se_beta", "est_se_beta",
                              "kappa", "mean_kappa", "se_kappa", "est_se_kappa")]

  library(xtable)
  write(
    print(xtable(to_print_ml), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("../../Draft/table/sim_negbin_sc", s, ".tex")
  )

  # re-order columns for table on moment estimate:
  to_print_moment <- to_print[, c("T",
                                  "tau", "mean_moment_tau", "se_moment_tau",
                                  "psi", "mean_moment_psi", "se_moment_psi",
                                  "beta", "mean_moment_beta", "se_moment_beta",
                                  "kappa", "mean_moment_kappa", "se_moment_kappa")]

  library(xtable)
  write(
    print(xtable(to_print_moment), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("../../Draft/table/sim_negbin_sc", s, "_moments.tex")
  )
}


