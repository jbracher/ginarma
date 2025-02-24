library(tidyverse)
library(xtable)

# True parameter values
vals_tau <- c(1, 1, 1)
vals_beta <- c(0.5, 0.2, 0.1)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_zeta <- c(0.5, 0.5, 0.5)

vals_lgt <- c(250, 500, 1000)  # Length of the trajectories
n_sim <- 1000  # Number of iterations

# Significance levels, for which the results shall be generated
significance <- c(0.05, 0.1)  

# Generate the table for INAR vs. INARMA ---------------------------------------
for (alpha in significance) {
  for (s in 1:3) {

    vals_lgt <- c(250, 500, 1000)

    scenario_temp <- matrix(nrow = 3, ncol = 2)
    rownames(scenario_temp) <- vals_lgt

    for (lgt in vals_lgt) {

      # Load the results
      res_null_vs_alternative <- read.csv(
        paste0(
          "LR_test/Results/INAR_fits_INARMA_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_null_vs_null <- read.csv(
        paste0(
          "LR_test/Results/INAR_fits_INAR_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_alternative_vs_null <- read.csv(
        paste0(
          "LR_test/Results/INARMA_fits_INAR_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_alternative_vs_alternative <- read.csv(
        paste0(
          "LR_test/Results/INARMA_fits_INARMA_holds_s", s, "_", lgt, ".csv"
        )
      )

      # Select only the convergent fits
      convergent_rejections <- res_null_vs_null$convergence == 0 &
        res_null_vs_alternative$convergence == 0
      convergent_power <- res_null_vs_null$convergence == 0 &
        res_alternative_vs_null$convergence == 0

      # Calculate the likelihood ratio statistic
      LR_stat <- -2 * (res_null_vs_null$max_lik[convergent_rejections] -
                         res_alternative_vs_null$max_lik[convergent_rejections])
      LR_stat_power <- -2 * (res_null_vs_alternative$max_lik[convergent_power] -
                               res_alternative_vs_alternative$max_lik[convergent_power])

      # Map the negative values of the LR-statistic to 0
      LR_stat_corrected <- ifelse(LR_stat > 0, LR_stat, 0)
      LR_stat_power_corrected <- ifelse(LR_stat_power > 0, LR_stat_power, 0)


        # Rejection rates

        # Compare with the quantile of the mixture of the chi^2(1) distribution
        # with point mass at 0
        scenario_temp[as.character(lgt), 1] <-
          sum(LR_stat_corrected > qchisq(1 - 2 * alpha, 1)) / sum(convergent_rejections)

        # Power

        # Compare with the quantile of the mixture of the chi^2(1) distribution
        # with point mass at 0
        scenario_temp[as.character(lgt), 2] <-
          sum(LR_stat_power_corrected > qchisq(1 - 2 * alpha, 1)) / sum(convergent_power)
    }
    # Format the results
    to_print <- cbind(s, vals_lgt, scenario_temp)
    to_print[, 2 + 1:ncol(scenario_temp)] <- format(round(100 * scenario_temp, digits = 2), nsmall = 2)
    to_print[grepl("NA", to_print)] <- "-"

    to_insert1 <- matrix(c("0.00", "", ""), nrow = 3, ncol = 1)
    to_insert2 <- matrix(c(format(vals_beta[s], nsmall = 2), "", ""), nrow = 3, ncol = 1)

    to_print <- cbind(to_print[, 1:2], to_insert1, to_print[, 3], to_insert2, to_print[, 4])

    write(
      print(
        xtable(to_print),
        only.contents = TRUE,
        include.rownames = FALSE, include.colnames = FALSE,
        hline.after = NULL),
      file = paste0("LR_test/Tables/LR_test_INAR_vs_INARMA_s", s, "_significance_", 100 * alpha, ".tex")
    )
  }
}

# Generate the table for INARMA vs. extended INARMA and INGARCH ----------------
for (alpha in significance) {
  for (s in 1:3) {

    vals_lgt <- c(250, 500, 1000)

    scenario_temp <- matrix(nrow = 3, ncol = 6)
    rownames(scenario_temp) <- vals_lgt

    for (lgt in vals_lgt) {

      # Load the results
      res_inarma_fits_vs_extended_holds <- read.csv(
        paste0(
          "LR_test/Results/INARMA_fits_extended_INARMA_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_inarma_fits_vs_inarma_holds <- read.csv(
        paste0(
          "LR_test/Results/INARMA_fits_INARMA_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_extended_fits_vs_inarma_holds <- read.csv(
        paste0(
          "LR_test/Results/extended_INARMA_fits_INARMA_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_extended_fits_vs_extended_holds <- read.csv(
        paste0(
          "LR_test/Results/extended_INARMA_fits_extended_INARMA_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_inarma_fits_vs_ingarch_holds <- read.csv(
        paste0(
          "LR_test/Results/INARMA_fits_INGARCH_holds_s", s, "_", lgt, ".csv"
        )
      )
      res_ingarch_fits_vs_ingarch_holds <- read.csv(
        paste0(
          "LR_test/Results/INGARCH_fits_INGARCH_holds_s", s, "_", lgt, ".csv"
        )
      )

      # Load the corrected thresholds (when available)
      corrected_thresholds <- read.csv(
        paste0(
          "LR_test/Results/INARMA_holds_correct_thresholds_s",
          s, "_", lgt, ".csv"
        )
      )

      corrected_thresholds_power <- read.csv(
        paste0(
          "LR_test/Results/extended_INARMA_holds_correct_thresholds_s",
          s, "_", lgt, ".csv"
        )
      )

      corrected_thresholds_power_vs_ingarch <- read.csv(
        paste0(
          "LR_test/Results/INGARCH_holds_correct_thresholds_s",
          s, "_", lgt, ".csv"
        )
      )

      # Select only the convergent fits
      convergent_rejections <- res_inarma_fits_vs_inarma_holds$convergence == 0 &
        res_extended_fits_vs_inarma_holds$convergence == 0 &
        !apply(res_inarma_fits_vs_inarma_holds[, 1:3], 1, anyNA) &
        !apply(res_extended_fits_vs_inarma_holds[, 1:3], 1, anyNA)
      convergent_power <- res_inarma_fits_vs_extended_holds$convergence == 0 &
        res_extended_fits_vs_extended_holds$convergence == 0 &
        !apply(res_inarma_fits_vs_extended_holds[, 1:3], 1, anyNA) &
        !apply(res_extended_fits_vs_extended_holds[, 1:3], 1, anyNA)
      convergent_power_vs_ingarch <- res_inarma_fits_vs_ingarch_holds$convergence == 0 &
        res_extended_fits_vs_extended_holds$convergence == 0 &
        !apply(res_inarma_fits_vs_ingarch_holds[, 1:3], 1, anyNA) &
        !apply(res_ingarch_fits_vs_ingarch_holds[, 1:3], 1, anyNA)

      # Calculate the likelihood ratio statistic
      LR_stat <- -2 * (res_inarma_fits_vs_inarma_holds$max_lik[convergent_rejections] -
                         res_extended_fits_vs_inarma_holds$max_lik[convergent_rejections])
      LR_stat_power <- -2 * (res_inarma_fits_vs_extended_holds$max_lik[convergent_power] -
                               res_extended_fits_vs_extended_holds$max_lik[convergent_power])
      LR_stat_power_vs_ingarch <- -2 * (res_inarma_fits_vs_ingarch_holds$max_lik[convergent_power_vs_ingarch] -
                                          res_ingarch_fits_vs_ingarch_holds$max_lik[convergent_power_vs_ingarch])

      # Map the negative values of the LR-statistic to 0
      LR_stat_corrected <- ifelse(LR_stat > 0, LR_stat, 0)
      LR_stat_power_corrected <- ifelse(LR_stat_power > 0, LR_stat_power, 0)
      LR_stat_power_vs_ingarch_corrected <- ifelse(LR_stat_power_vs_ingarch > 0,
                                                   LR_stat_power_vs_ingarch, 0)

      # Rejection rates

      # Compare with the quantile of the mixture of the chi^2(1) distribution
      # with point mass at 0
      scenario_temp[as.character(lgt), 1] <-
        sum(LR_stat_corrected > qchisq(1 - 2 * alpha, 1)) / sum(convergent_rejections)
      # Compare with the corrected thresholds
      scenario_temp[as.character(lgt), 2] <-
        sum(LR_stat_corrected > corrected_thresholds[convergent_rejections, paste0("critical_", alpha)]) / sum(convergent_rejections)

      # Power (INARMA vs. extended INARMA)

      # Compare with the quantile of the mixture of the chi^2(1) distribution
      # with point mass at 0
      scenario_temp[as.character(lgt), 3] <-
        sum(LR_stat_power_corrected > qchisq(1 - 2 * alpha, 1)) / sum(convergent_power)
      # Compare with the corrected thresholds
      scenario_temp[as.character(lgt), 4] <-
        sum(LR_stat_power_corrected > corrected_thresholds_power[convergent_power, paste0("critical_", alpha)]) / sum(convergent_power)

      # Power (INARMA vs. INGARCH)

      # Compare with the quantile of the mixture of the chi^2(1) distribution
      # with point mass at 0
      scenario_temp[as.character(lgt), 5] <-
        sum(LR_stat_power_vs_ingarch_corrected > qchisq(1 - 2 * alpha, 1)) / sum(convergent_power)
      # Compare with the corrected thresholds
      scenario_temp[as.character(lgt), 6] <-
        sum(LR_stat_power_vs_ingarch_corrected > corrected_thresholds_power_vs_ingarch[convergent_power_vs_ingarch, paste0("critical_", alpha)]) / sum(convergent_power)
    }
    # Format the results
    to_print <- cbind(s, vals_lgt, scenario_temp)
    to_print[, 2 + 1:ncol(scenario_temp)] <- format(round(100 * scenario_temp, digits = 2), nsmall = 2, width = 6)
    to_print[grepl("NA", to_print)] <- "-"

    to_insert1 <- matrix(c("0.00", "", ""), nrow = 3, ncol = 1)
    to_insert2 <- matrix(c("0.50", "", ""), nrow = 3, ncol = 1)
    to_insert3 <- matrix(c("1.00", "", ""), nrow = 3, ncol = 1)

    to_print <- cbind(to_print[, 1:2], to_insert1, to_print[, 3:4], to_insert2,
                      to_print[, 5:6], to_insert3, to_print[, 7:8])

    write(
      print(
        xtable(to_print, ),
        only.contents = TRUE,
        include.rownames = FALSE, include.colnames = FALSE,
        hline.after = NULL),
      file = paste0("LR_test/Tables/LR_test_INARMA_vs_extended_INARMA_and_INGARCH_s", s, "_significance_", 100 * alpha, ".tex")
    )
  }
}
