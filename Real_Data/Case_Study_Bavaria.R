# Case study on measles and mumps in Bavaria

library(inarma)

# load data:
data("measles")
data("mumps")
list_data <- list(mumps, measles)
names_diseases <- c("measles", "mumps")
names_units_to_print <- c("Measles", "Mumps")

# descriptive plot of data:
pdf("Real_Data/Plots/data.pdf", width = 8, height = 4)
layout_matr <- matrix(c(1, 1, 4, 4,
                        2, 3, 5, 6), byrow = TRUE, nrow = 2)

layout(layout_matr)

ylims_ts <- list(measles = c(0, 15), mumps = c(0, 15))
ylims_tab <- list(measles = c(0, 120), mumps = c(0, 120))
xlims_tab <- list(measles = c(0, 15), mumps = c(0, 15))

for(i in seq_along(names_diseases)){
  disease <- names_diseases[i]
  data_temp <- get(disease)
  data_temp$time <- data_temp$year + data_temp$week/52

  par(mar = c(1, 4, 3, 1), las = 1)
  plot(data_temp$time, data_temp$value, type = "l", xlab = "", ylab = "cases",
       main = names_units_to_print[i],
       axes = FALSE, ylim = ylims_ts[[disease]])
  axis(1, at = 2014:2019 + 0.5, labels = 2014:2019, tick = FALSE)
  axis(1, labels = FALSE)
  axis(2)
  box()

  par(mar = c(4, 4, 3, 1), las = 1)
  plot(table(data_temp$value), xlab = "cases", ylab = "frequency",
       xlim = xlims_tab[[disease]], ylim = ylims_tab[[disease]], axes = FALSE)
  axis(2)
  axis(1, at = c(0, 5, 10, 15))

  acf(data_temp$value, lag.max = 10, ci.col = "black", main = "", axes = FALSE,
      xlab = "lag")
  axis(1)
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1))
  box()
  text(4, 0.84, paste("mean =", round(mean(data_temp$value), 2)))
  text(4, 0.61, paste("   var =", round(var(data_temp$value), 2)))
}
dev.off()

#############
## fit models
#############

### Measles
fits_measles <- list()

# fit inarch models
fits_measles$inarch <- fit_inarch(measles$value, family = "Poisson")
fits_measles$hinarch <- fit_inarch(measles$value, family = "Hermite")
fits_measles$nbinarch <- fit_inarch(measles$value, family = "NegBin")

# fit ingarch models
fits_measles$ingarch <- fit_ingarch(measles$value, family = "Poisson")
fits_measles$hingarch <- fit_ingarch(measles$value, family = "Hermite", control_optim = list(maxit = 800))
fits_measles$nbingarch <- fit_ingarch(measles$value, family = "NegBin", control_optim = list(maxit = 800))

# fit inar models
fits_measles$inar <- fit_inar(measles$value, family = "Poisson", offspring = "binomial")
fits_measles$hinar <- fit_inar(measles$value, family = "Hermite", offspring = "binomial")
fits_measles$nbinar <- fit_inar(measles$value, family = "NegBin", offspring = "binomial")

# fit inarma models (using Poisson initialization):
fits_measles$inarma <- fit_inarma(measles$value, family = "Poisson", offspring = "binomial")
fits_measles$hinarma <- fit_inarma(measles$value, family = "Hermite", offspring = "binomial", control_optim = list(maxit = 800))
fits_measles$nbinarma <- fit_inarma(measles$value, family = "NegBin", offspring = "binomial", control_optim = list(maxit = 1000))

# print AIC values:
old_measles_AIC <- c(
  inarch = 1159.13, hinarch = 1082.40, nbinarch = 1055.04,
  ingarch = 1096.91, hingarch = 1046.09, nbingarch = 1028.23,
  inar = 1232.94, hinar = 1122.68, nbinar = 1068.77,
  inarma = 1166.26, hinarma = 1094.07, nbinarma = 1046.65
)
unlist(lapply(fits_measles, FUN = function(mod) mod$AIC))
unlist(lapply(fits_measles, FUN = function(mod) mod$opt$convergence))


### Mumps
fits_mumps <- list()

# fit inarch models
fits_mumps$inarch <- fit_inarch(mumps$value, family = "Poisson")
fits_mumps$hinarch <- fit_inarch(mumps$value, family = "Hermite")
fits_mumps$nbinarch <- fit_inarch(mumps$value, family = "NegBin")

# fit ingarch models
fits_mumps$ingarch <- fit_ingarch(mumps$value, family = "Poisson")
fits_mumps$hingarch <- fit_ingarch(mumps$value, family = "Hermite", control_optim = list(maxit = 800))
fits_mumps$nbingarch <- fit_ingarch(mumps$value, family = "NegBin")

# fit inar models (using Poisson initialization):
fits_mumps$inar <- fit_inar(mumps$value, family = "Poisson", offspring = "binomial")
fits_mumps$hinar <- fit_inar(mumps$value, family = "Hermite", offspring = "binomial")
fits_mumps$nbinar <- fit_inar(mumps$value, family = "NegBin", offspring = "binomial")

# fit inarma models (using Poisson initialization):
fits_mumps$inarma <- fit_inarma(mumps$value, family = "Poisson", offspring = "binomial")
fits_mumps$hinarma <- fit_inarma(mumps$value, family = "Hermite", offspring = "binomial", control_optim = list(maxit = 800))
fits_mumps$nbinarma <- fit_inarma(mumps$value, family = "NegBin", offspring = "binomial", control_optim = list(maxit = 800))

# print AIC values:
old_mumps_AIC <- c(
  inarch = 1274.26, hinarch = 1249.33, nbinarch = 1244.75,
  ingarch = 1238.27, hingarch = 1224.43, nbingarch = 1222.86,
  inar = 1283.22, hinar = 1252.77, nbinar = 1245.57,
  inarma = 1257.34, hinarma = 1235.48, nbinarma = 1231.73
)
unlist(lapply(fits_mumps, FUN = function(mod) mod$AIC))
unlist(lapply(fits_mumps, FUN =  function(mod) mod$convergence))

# print AIC values:
old_mumps_AIC <- c(
  inarch = 1274.26, hinarch = 1249.33, nbinarch = 1244.75,
  ingarch = 1238.27, hingarch = 1224.43, nbingarch = 1222.86,
  inar = 1283.22, hinar = 1252.77, nbinar = 1245.57,
  inarma = 1257.34, hinarma = 1235.48, nbinarma = 1231.73
)

# save results:
save(fits_measles, fits_mumps, file = "Real_Data/Fits/fits_measles_mumps_bavaria.rda")

#######################
## create output tables
#######################

### Parameter scale ------------------------------------------------------------

# INAR

# initialize table:
tab_inar <- data.frame(
  model = c("inar", "hinar", "nbinar"),
  nu_measles = NA,
  alpha_measles = NA,
  beta_measles = NA,
  psi_measles = NA,
  BIC_measles = NA,
  nu_mumps = NA,
  alpha_mumps = NA,
  beta_mumps = NA,
  psi_mumps = NA,
  BIC_mumps = NA
)

# fill table:
for(i in 1:3){
  mod <- tab_inar$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inar[i, "nu_measles"] <- fit_temp$coefficients["tau"]
  tab_inar[i, "alpha_measles"] <- fit_temp$coefficients["kappa"]
  tab_inar[i, "psi_measles"] <- fit_temp$coefficients["psi"]
  tab_inar[i, "BIC_measles"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood

  fit_temp <- fits_mumps[[mod]]
  tab_inar[i, "nu_mumps"] <- fit_temp$coefficients["tau"]
  tab_inar[i, "alpha_mumps"] <- fit_temp$coefficients["kappa"]
  tab_inar[i, "psi_mumps"] <- fit_temp$coefficients["psi"]
  tab_inar[i, "BIC_mumps"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood
}
tab_inar$model <- toupper(tab_inar$model)
tab_inar_text <- data.frame(sapply(
  tab_inar,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_inar_text$psi_measles[1] <- "-"
tab_inar_text$psi_mumps[1] <- "-"
tab_inar_text$beta_measles[1:3] <- "-"
tab_inar_text$beta_mumps[1:3] <- "-"

# INARCH

# initialize table:
tab_inarch <- data.frame(
  model = c("inarch", "hinarch", "nbinarch"),
  nu_measles = NA,
  alpha_measles = NA,
  beta_measles = NA,
  psi_measles = NA,
  BIC_measles = NA,
  nu_mumps = NA,
  alpha_mumps = NA,
  beta_mumps = NA,
  psi_mumps = NA,
  BIC_mumps = NA
)
# fill table:
for(i in 1:3){
  mod <- tab_inarch$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inarch[i, "nu_measles"] <- fit_temp$coefficients["tau"] *
    ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_inarch[i, "alpha_measles"] <- fit_temp$coefficients["kappa"] *
    ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_inarch[i, "psi_measles"] <- fit_temp$psi["psi"]
  tab_inarch[i, "BIC_measles"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood

  fit_temp <- fits_mumps[[mod]]
  tab_inarch[i, "nu_mumps"] <- fit_temp$coefficients["tau"] *
    ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_inarch[i, "alpha_mumps"] <- fit_temp$coefficients["kappa"] *
    ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_inarch[i, "psi_mumps"] <- fit_temp$psi["psi"]
  tab_inarch[i, "BIC_mumps"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood
}
tab_inarch$model <- c("Poisson INARCH", "Hermite INARCH", "NegBin INARCH")
tab_inarch_text <- data.frame(sapply(
  tab_inarch,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_inarch_text$psi_measles[1] <- "-"
tab_inarch_text$psi_mumps[1] <- "-"
tab_inarch_text$beta_measles[1:3] <- "-"
tab_inarch_text$beta_mumps[1:3] <- "-"


# INARMA

# initialize table:
tab_inarma <- data.frame(
  model = c("inarma", "hinarma", "nbinarma"),
  nu_measles = NA,
  kappa_measles = NA,
  beta_measles = NA,
  psi_measles = NA,
  BIC_measles = NA,
  nu_mumps = NA,
  kappa_mumps = NA,
  beta_mumps = NA,
  psi_mumps = NA,
  BIC_mumps = NA
)
# fill table:
for(i in 1:3){
  mod <- tab_inarma$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inarma[i, "nu_measles"] <- fit_temp$coefficients["tau"]
  tab_inarma[i, "kappa_measles"] <- fit_temp$coefficients["kappa"]
  tab_inarma[i, "beta_measles"] <- fit_temp$coefficients["beta"]
  tab_inarma[i, "psi_measles"] <- fit_temp$coefficients["psi"]
  tab_inarma[i, "BIC_measles"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood

  fit_temp <- fits_mumps[[mod]]
  tab_inarma[i, "nu_mumps"] <- fit_temp$coefficients["tau"]
  tab_inarma[i, "kappa_mumps"] <- fit_temp$coefficients["kappa"]
  tab_inarma[i, "beta_mumps"] <- fit_temp$coefficients["beta"]
  tab_inarma[i, "psi_mumps"] <- fit_temp$coefficients["psi"]
  tab_inarma[i, "BIC_mumps"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood
}
tab_inarma$model <- c("Poisson INARMA", "Hermite INARMA", "NegBin INARMA")
tab_inarma_text <- data.frame(sapply(
  tab_inarma,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_inarma_text$psi_measles[1] <- "-"
tab_inarma_text$psi_mumps[1] <- "-"

# INGARCH

# initialize table:
tab_ingarch <- data.frame(
  model = c("ingarch", "hingarch", "nbingarch"),
  nu_measles = NA,
  alpha_measles = NA,
  beta_measles = NA,
  psi_measles = NA,
  BIC_measles = NA,
  nu_mumps = NA,
  alpha_mumps = NA,
  beta_mumps = NA,
  psi_mumps = NA,
  BIC_mumps = NA
)
# fill table:
for(i in 1:3){
  mod <- tab_ingarch$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_ingarch[i, "nu_measles"] <- fit_temp$coefficients["tau"] *
    (1 - fit_temp$coefficients["beta"]) * ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_ingarch[i, "alpha_measles"] <- fit_temp$coefficients["kappa"] *
    (1 - fit_temp$coefficients["beta"]) * ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_ingarch[i, "beta_measles"] <- fit_temp$coefficients["beta"]
  tab_ingarch[i, "psi_measles"] <- fit_temp$psi["psi"]
  tab_ingarch[i, "BIC_measles"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood

  fit_temp <- fits_mumps[[mod]]
  tab_ingarch[i, "nu_mumps"] <- fit_temp$coefficients["tau"] *
    (1 - fit_temp$coefficients["beta"]) * ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_ingarch[i, "alpha_mumps"] <- fit_temp$coefficients["kappa"] *
    (1 - fit_temp$coefficients["beta"]) * ifelse(i == 1, 1, fit_temp$coefficients["theta"])
  tab_ingarch[i, "beta_mumps"] <- fit_temp$coefficients["beta"]
  tab_ingarch[i, "psi_mumps"] <- fit_temp$psi["psi"]
  tab_ingarch[i, "BIC_mumps"] <- fit_temp$dim * log(fit_temp$nobs) -
    2 * fit_temp$loglikelihood
}
tab_ingarch$model <- c("Poisson INGARCH", "Hermite INGARCH", "NegBin INGARCH")
tab_ingarch_text <- data.frame(sapply(
  tab_ingarch,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_ingarch_text$psi_measles[1] <- "-"
tab_ingarch_text$psi_mumps[1] <- "-"

# write out:
library(xtable)

tab_inar <- tab_inar_text
tab_inarch <- tab_inarch_text
tab_ingarch <- tab_ingarch_text
tab_inarma <- tab_inarma_text

for(tab in c("tab_inar", "tab_inarch", "tab_inarma", "tab_ingarch")){
  tab_temp <- get(tab)

  write(
    print(xtable(tab_temp), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("Real_Data/Tables/", tab, ".tex")
  )
}


# Epidemiological scales -------------------------------------------------------


# INAR

# initialize table:
tab_inar_epi <- tab_se_inar_epi <- data.frame(model = c("inar", "hinar", "nbinar"),
                                              immig_measles = NA, R_measles = NA, GT_measles = NA, CS_measles = NA, AIC_measles = NA,
                                              immig_mumps = NA, R_mumps = NA,  GT_mumps = NA, CS_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inar_epi$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inar_epi[i, "R_measles"] <- fit_temp$coefficients["kappa"]
  tab_se_inar_epi[i, "R_measles"] <- fit_temp$se["kappa"]

  tab_inar_epi[i, "immig_measles"] <- fit_temp$coefficients["tau"]
  tab_se_inar_epi[i, "immig_measles"] <- fit_temp$se["tau"]

  # generation time is 1 by construction:
  tab_inar_epi[i, "GT_measles"] <- 1
  # cluster size is 1 by construction:
  tab_inar_epi[i, "CS_measles"] <- 1
  tab_inar_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inar_epi[i, "R_mumps"] <- fit_temp$coefficients["kappa"]
  tab_se_inar_epi[i, "R_mumps"] <- fit_temp$se["kappa"]

  tab_inar_epi[i, "immig_mumps"] <- fit_temp$coefficients["tau"]
  tab_se_inar_epi[i, "immig_mumps"] <- fit_temp$se["tau"]
  # generation time is 1 by construction:
  tab_inar_epi[i, "GT_mumps"] <- 1
  # cluster size is 1 by construction:
  tab_inar_epi[i, "CS_mumps"] <- 1
  tab_inar_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inar_epi$model <- c("Poisson INAR", "Hermite INAR", "NegBin INAR")
tab_inar_epi_text <- data.frame(sapply(
  tab_inar_epi,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_inar_epi_text$CS_measles[1:3] <- "1*"
tab_inar_epi_text$CS_mumps[1:3] <- "1*"
tab_inar_epi_text$GT_measles[1:3] <- "1*"
tab_inar_epi_text$GT_mumps[1:3] <- "1*"

# INARCH

# initialize table:
tab_inarch_epi <- data.frame(model = c("inarch", "hinarch", "nbinarch"),
                             immig_measles = NA, R_measles = NA, GT_measles = NA, CS_measles = NA, AIC_measles = NA,
                             immig_mumps = NA, R_mumps = NA, GT_mumps = NA, CS_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inarch_epi$model[i]
  family <- switch(mod, inarch = "Poisson", hinarch = "Hermite", nbinarch = "NegBin")

  fit_temp <- fits_measles[[mod]]
  tab_inarch_epi[i, "R_measles"] <- fit_temp$coefficients["kappa"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  tab_inarch_epi[i, "immig_measles"] <- fit_temp$coefficients["tau"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  # generation time is 1 by construction:
  tab_inarch_epi[i, "GT_measles"] <- 1
  tab_inarch_epi[i, "CS_measles"] <- fit_temp$coefficients["theta"]
  tab_inarch_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inarch_epi[i, "R_mumps"] <- fit_temp$coefficients["kappa"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  tab_inarch_epi[i, "immig_mumps"] <- fit_temp$coefficients["tau"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  # generation time is 1 by construction:
  tab_inarch_epi[i, "GT_mumps"] <- 1
  tab_inarch_epi[i, "CS_mumps"] <- fit_temp$coefficients["theta"]
  tab_inarch_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inarch_epi$model <- c("Poisson INARCH", "Hermite INARCH", "NegBin INARCH")
tab_inarch_epi_text <- data.frame(sapply(
  tab_inarch_epi,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_inarch_epi_text$CS_measles[1] <- "1*"
tab_inarch_epi_text$CS_mumps[1] <- "1*"
tab_inarch_epi_text$GT_measles[1:3] <- "1*"
tab_inarch_epi_text$GT_mumps[1:3] <- "1*"

# INARMA

# initialize table:
tab_inarma_epi <- data.frame(model = c("inarma", "hinarma", "nbinarma"),
                             immig_measles = NA, R_measles = NA, GT_measles = NA, CS_measles = NA, AIC_measles = NA,
                             immig_mumps = NA, R_mumps = NA, GT_mumps = NA, CS_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inarma_epi$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inarma_epi[i, "immig_measles"] <- fit_temp$coefficients["tau"]
  tab_inarma_epi[i, "R_measles"] <- fit_temp$coefficients["kappa"]
  tab_inarma_epi[i, "GT_measles"] <- 1 / (1 - fit_temp$coefficients["beta"])
  # cluster size is 1 by construction
  tab_inarma_epi[i, "CS_measles"] <- 1
  tab_inarma_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inarma_epi[i, "immig_mumps"] <- fit_temp$coefficients["tau"]
  tab_inarma_epi[i, "R_mumps"] <- fit_temp$coefficients["kappa"]
  tab_inarma_epi[i, "GT_mumps"] <- 1 / (1 - fit_temp$coefficients["beta"])
  # cluster size is 1 by construction
  tab_inarma_epi[i, "CS_mumps"] <- 1
  tab_inarma_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inarma_epi$model <- c("Poisson INARMA", "Hermite INARMA", "NegBin INARMA")
tab_inarma_epi_text <- data.frame(sapply(
  tab_inarma_epi,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
))
tab_inarma_epi_text$CS_measles[1:3] <- "1*"
tab_inarma_epi_text$CS_mumps[1:3] <- "1*"


# INGARCH:
# initialize table:
tab_ingarch_epi <- data.frame(model = c("ingarch", "hingarch", "nbingarch"),
                              immig_measles = NA, R_measles = NA, GT_measles = NA, CS_measles = NA, AIC_measles = NA,
                              immig_mumps = NA, R_mumps = NA, GT_mumps = NA, CS_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_ingarch_epi$model[i]
  family <- switch(mod, ingarch = "Poisson", hingarch = "Hermite", nbingarch = "NegBin")

  fit_temp <- fits_measles[[mod]]
  tab_ingarch_epi[i, "R_measles"] <- fit_temp$coefficients["kappa"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  tab_ingarch_epi[i, "immig_measles"] <- fit_temp$coefficients["tau"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  tab_ingarch_epi[i, "GT_measles"] <- 1/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "CS_measles"] <- fit_temp$coefficients["theta"]
  tab_ingarch_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_ingarch_epi[i, "R_mumps"] <- fit_temp$coefficients["kappa"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  tab_ingarch_epi[i, "immig_mumps"] <- fit_temp$coefficients["tau"] *
    ifelse(family == "Poisson", 1, fit_temp$coefficients["theta"])
  tab_ingarch_epi[i, "GT_mumps"] <- 1/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "CS_mumps"] <- fit_temp$coefficients["theta"]
  tab_ingarch_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_ingarch_epi$model <- c("Poisson INGARCH", "Hermite INGARCH", "NegBin INGARCH")

tab_ingarch_epi_text <- as.data.frame(sapply(
  tab_ingarch_epi,
  function (x) {as.character(format(x, digits = 2, nsmall = 2, scientific = FALSE))}
  ))
tab_ingarch_epi_text$CS_measles[1] <- "1*"
tab_ingarch_epi_text$CS_mumps[1] <- "1*"



# write out:
library(xtable)

tab_inar_epi <- tab_inar_epi_text
tab_inarch_epi <- tab_inarch_epi_text
tab_ingarch_epi <- tab_ingarch_epi_text
tab_inarma_epi <- tab_inarma_epi_text

for(tab in c("tab_inar_epi", "tab_inarch_epi", "tab_inarma_epi", "tab_ingarch_epi")){
  tab_temp <- get(tab)

  write(
    print(xtable(tab_temp), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("Real_Data/Tables/", tab, ".tex")
  )
}

#####################################
###### diagnostic plot: residual ACF

# define colours:
library(RColorBrewer)

# get fits (in case this code is run without re-doing the fit)
# load("R/Real_Data/Fits/fits_measles_mumps_bavaria.rda")

# compute residual autocorrelations and variances:
acs_measles <- lapply(fits_measles, function(x) acf(x$pearson_residuals))
vars_measles <- sapply(fits_measles, function(x) var(x$pearson_residuals))
acs_mumps <- lapply(fits_mumps, function(x) acf(x$pearson_residuals))
vars_mumps <- sapply(fits_mumps, function(x) var(x$pearson_residuals))

# plot for INARMA:
pdf("Real_Data/Plots/residuals_inarma.pdf", width = 8, height = 3.5)

# define colours
cols <- brewer.pal(6, "Paired")
names(cols) <- c("inar", "inarma",
                 "hinar", "hinarma",
                 "nbinar", "nbinarma")

# order in which results are plotted
order_models <- c("inar", "inarma",
                  "hinar", "hinarma",
                  "nbinar", "nbinarma")

# structure plot area
layout_matr <- matrix(c(1, 2, 2, 2, 5, 3, 4, 4, 4, 5), nrow = 2, byrow = TRUE)
layout(layout_matr)
par(las = 1, mar = c(4, 5, 3, 1))

# measles:
plot(rep(0:1, 3), vars_measles[order_models], xlim = c(-0.5, 1.5),
     pch = 15, col = cols[order_models], ylim = c(0, 3), axes = FALSE,
     xlab = "", ylab = "variance of \n Pearson residuals")
axis(2)
box()

vals_x <- 1:6
plot(vals_x - 0.25, acs_measles$inar$acf[2:7], type = "h", ylim = c(0, 0.4), xlim = c(0.5, 6.5),
     col = cols["inar"], lwd = 2, xlab = "lag", ylab = "residual autocorr.",
     main = "Measles")
abline(h = 2/sqrt(length(fits_measles$inarch$fitted_values)), lty = "dashed")

points(vals_x - 0.15, acs_measles$hinar$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["hinar"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x - 0.05, acs_measles$nbinar$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["nbinar"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.05, acs_measles$inarma$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["inarma"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.15, acs_measles$hinarma$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["hinarma"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.25, acs_measles$nbinarma$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["nbinarma"], lwd = 2, xlab = "lag", ylab = "ACF")

# mumps:
plot(rep(0:1, 3), vars_mumps[order_models], xlim = c(-0.5, 1.5),
     pch = 15, col = cols[order_models], ylim = c(0, 3), axes = FALSE,
     xlab = "", ylab = "variance of \n Pearson residuals")
axis(2)
box()

vals_x <- 1:6
plot(vals_x - 0.25, acs_mumps$inar$acf[2:7], type = "h", ylim = c(0, 0.4), xlim = c(0.5, 6.5),
     col = cols["inar"], lwd = 2, xlab = "lag", ylab = "residual autocorr.",
     main = "Mumps")
abline(h = 2/sqrt(length(fits_mumps$inarch$fitted_values)), lty = "dashed")

points(vals_x - 0.15, acs_mumps$hinar$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["hinar"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x - 0.05, acs_mumps$nbinar$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["nbinar"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.05, acs_mumps$inarma$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["inarma"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.15, acs_mumps$hinarma$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["hinarma"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.25, acs_mumps$nbinarma$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols["nbinarma"], lwd = 2, xlab = "lag", ylab = "ACF")

par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab ="", ylab = "")
order_models_legend <- c("inar", "hinar", "nbinar",
                         "inarma", "hinarma", "nbinarma")

legend("right", col = cols[order_models_legend],
       legend = c("Poisson INAR", "Hermite INAR", "NegBin INAR",
                  "Poisson INARMA", "Hermite INARMA", "NegBin INARMA"), lwd = 2, bty = "n")
dev.off()



# plot for INGARCH:
pdf("Real_Data/Plots/residuals_ingarch.pdf", width = 8, height = 3.5)

# define colours
cols_inarch <- brewer.pal(6, "Paired")
names(cols_inarch) <- c("inarch", "ingarch",
                        "hinarch", "hingarch",
                        "nbinarch", "nbingarch")

# order in which results are plotted
order_models <- c("inarch", "ingarch",
                  "hinarch", "hingarch",
                  "nbinarch", "nbingarch")

# structure plot area
layout_matr <- matrix(c(1, 2, 2, 2, 5, 3, 4, 4, 4, 5), nrow = 2, byrow = TRUE)
layout(layout_matr)
par(las = 1, mar = c(4, 5, 3, 1))

# measles:
plot(rep(0:1, 3), vars_measles[order_models], xlim = c(-0.5, 1.5),
     pch = 15, col = cols_inarch[order_models], ylim = c(0, 3), axes = FALSE,
     xlab = "", ylab = "variance of \n Pearson residuals")
axis(2)
box()

vals_x <- 1:6
plot(vals_x - 0.25, acs_measles$inarch$acf[2:7], type = "h", ylim = c(0, 0.4), xlim = c(0.5, 6.5),
     col = cols_inarch["inarch"], lwd = 2, xlab = "lag", ylab = "residual autocorr.",
     main = "Measles")
abline(h = 2/sqrt(length(fits_measles$inarch$fitted_values)), lty = "dashed")

points(vals_x - 0.15, acs_measles$hinarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["hinarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x - 0.05, acs_measles$nbinarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["nbinarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.05, acs_measles$ingarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["ingarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.15, acs_measles$hingarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["hingarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.25, acs_measles$nbingarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["nbingarch"], lwd = 2, xlab = "lag", ylab = "ACF")

# mumps:

plot(rep(0:1, 3), vars_mumps[order_models], xlim = c(-0.5, 1.5),
     pch = 15, col = cols_inarch[order_models], ylim = c(0, 3), axes = FALSE,
     xlab = "", ylab = "variance of \n Pearson residuals")
axis(2)
box()

vals_x <- 1:6
plot(vals_x - 0.25, acs_mumps$inarch$acf[2:7], type = "h", ylim = c(0, 0.4), xlim = c(0.5, 6.5),
     col = cols_inarch["inarch"], lwd = 2, xlab = "lag", ylab = "residual autocorr.",
     main = "Mumps")
abline(h = 2/sqrt(length(fits_mumps$inarch$fitted_values)), lty = "dashed")

points(vals_x - 0.15, acs_mumps$hinarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["hinarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x - 0.05, acs_mumps$nbinarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["nbinarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.05, acs_mumps$ingarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["ingarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.15, acs_mumps$hingarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["hingarch"], lwd = 2, xlab = "lag", ylab = "ACF")

points(vals_x + 0.25, acs_mumps$nbingarch$acf[2:7], type = "h", ylim = c(0, 0.5),
       col = cols_inarch["nbingarch"], lwd = 2, xlab = "lag", ylab = "ACF")

par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab ="", ylab = "")
order_models_legend <- c("inarch", "hinarch", "nbinarch",
                         "ingarch", "hingarch", "nbingarch")

legend("right", col = cols_inarch[order_models_legend],
       legend = c("Poisson INARCH", "Hermite INARCH", "NegBin INARCH",
                  "Poisson INGARCH", "Hermite INGARCH", "NegBin INGARCH"), lwd = 2, bty = "n")
dev.off()


#####################################
###### diagnostic plot: fits

pdf("Real_Data/Plots/fits_measles_inar.pdf", width = 8, height = 4.5)
par(mar = c(4, 4, 2, 1), mfrow = c(2, 3), las = 1)
plot(fits_measles$inar, type = "fit", interval_level = 0.9, legend = FALSE)
title("Poisson INAR", line = 1)
plot(fits_measles$hinar, type = "fit", interval_level = 0.9, legend = FALSE)
title("Hermite INAR", line = 1)
plot(fits_measles$nbinar, type = "fit", interval_level = 0.9, legend = FALSE)
title("Negative binomial INAR", line = 1)

plot(fits_measles$inarma, type = "fit", interval_level = 0.9, legend = FALSE)
title("Poisson INARMA", line = 1)
plot(fits_measles$hinarma, type = "fit", interval_level = 0.9, legend = FALSE)
title("Hermite INARMA", line = 1)
plot(fits_measles$nbinarma, type = "fit", interval_level = 0.9, legend = FALSE)
title("Negative binomial INARMA", line = 1)
dev.off()



pdf("Real_Data/Plots/fits_mumps_inar.pdf", width = 8, height = 4.5)
par(mar = c(4, 4, 2, 1), mfrow = c(2, 3), las = 1)
plot(fits_mumps$inar, type = "fit", interval_level = 0.9, legend = FALSE)
title("Poisson INAR", line = 1)
plot(fits_mumps$hinar, type = "fit", interval_level = 0.9, legend = FALSE)
title("Hermite INAR", line = 1)
plot(fits_mumps$nbinar, type = "fit", interval_level = 0.9, legend = FALSE)
title("Negative binomial INAR", line = 1)

plot(fits_mumps$inarma, type = "fit", interval_level = 0.9, legend = FALSE)
title("Poisson INARMA", line = 1)
plot(fits_mumps$hinarma, type = "fit", interval_level = 0.9, legend = FALSE)
title("Hermite INARMA", line = 1)
plot(fits_mumps$nbinarma, type = "fit", interval_level = 0.9, legend = FALSE)
title("Negative binomial INARMA", line = 1)
dev.off()
