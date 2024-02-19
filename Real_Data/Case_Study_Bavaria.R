####################################################
# Case study on measles and mumps in Bavaria
# Johannes Bracher, johannes.bracher@kit.edu
####################################################

#####################################
###### preparations

# load package
library(inarma)

# set working directory to current folder
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get some additional functions to fit INGARCH models etc
source("helper_functions.R")

# load data:
data("measles")
data("mumps")
list_data <- list(mumps, measles)
names_diseases <- c("measles", "mumps")
names_units_to_print <- c("Measles", "Mumps")

# descriptive plot of data:
pdf("Plots/data.pdf", width = 8, height = 4)
layout_matr <- matrix(c(1, 1, 4, 4,
                        2, 3, 5, 6), byrow = TRUE, nrow = 2)

layout(layout_matr)
par(mar = c(2, 4, 3, 1), las = 1)

ylims_ts <- list(measles = c(0, 15), mumps = c(0, 15))
ylims_tab <- list(measles = c(0, 120), mumps = c(0, 120))
xlims_tab <- list(measles = c(0, 15), mumps = c(0, 15))



for(i in seq_along(names_diseases)){
  disease <- names_diseases[i]
  data_temp <- get(disease)
  data_temp$time <- data_temp$year + data_temp$week/52
  plot(data_temp$time, data_temp$value, type = "l", xlab = "", ylab = "cases",
       main = names_units_to_print[i],
       axes = FALSE, ylim = ylims_ts[[disease]])
  axis(1, at = 2014:2019 + 0.5, labels = 2014:2019, tick = FALSE)
  axis(1, labels = FALSE)
  axis(2)
  box()

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

#####################################
###### fit models

### measles
fits_measles <- list()

# fit inarch models
fits_measles$inarch <- fit_inarch(measles$value)
fits_measles$hinarch <- fit_hinarch(measles$value)
fits_measles$nbinarch <- fit_nbinarch2(measles$value)

# fit ingarch models
fits_measles$ingarch <- fit_ingarch(measles$value)
fits_measles$hingarch <- fit_hingarch(measles$value)
fits_measles$nbingarch <- fit_nbingarch2(measles$value)

# fit inar models (using Poisson initialization):
fits_measles$inar <- fit_inar(measles$value, family = "Poisson")
fits_measles$hinar <- fit_inar(measles$value, family = "Hermite")
fits_measles$nbinar <- fit_inar(measles$value, family = "NegBin")

# fit inarma models (using Poisson initialization):
fits_measles$inarma <- fit_inarma(measles$value, family = "Poisson")
fits_measles$hinarma <- fit_inarma(measles$value, family = "Hermite")
fits_measles$nbinarma <- fit_inarma(measles$value, family = "NegBin")

# print AIC values:
lapply(fits_measles, FUN = function(mod) mod$AIC)
lapply(fits_measles, FUN = function(mod) mod$opt$convergence)


### mumps
fits_mumps <- list()

# fit inarch models
fits_mumps$inarch <- fit_inarch(mumps$value)
fits_mumps$hinarch <- fit_hinarch(mumps$value)
fits_mumps$nbinarch <- fit_nbinarch2(mumps$value)

# fit ingarch models
fits_mumps$ingarch <- fit_ingarch(mumps$value)
fits_mumps$hingarch <- fit_hingarch(mumps$value)
fits_mumps$nbingarch <- fit_nbingarch2(mumps$value)

# fit inar models (using Poisson initialization):
fits_mumps$inar <- fit_inar(mumps$value, family = "Poisson")
fits_mumps$hinar <- fit_inar(mumps$value, family = "Hermite")
fits_mumps$nbinar <- fit_inar(mumps$value, family = "NegBin")

# fit inarma models (using Poisson initialization):
fits_mumps$inarma <- fit_inarma(mumps$value, family = "Poisson")
fits_mumps$hinarma <- fit_inarma(mumps$value, family = "Hermite")
fits_mumps$nbinarma <- fit_inarma(mumps$value, family = "NegBin")

# print AIC values:
lapply(fits_mumps, FUN =  function(mod) mod$AIC)

# save results:
save(fits_measles, fits_mumps, file = "Fits/fits_measles_mumps_bavaria.rda")

load("Fits/fits_measles_mumps_bavaria.rda")

#####################################
###### create output table:

#####################################
###### regular parameters:

# INAR:
# initialize table:
tab_inar <- data.frame(model = c("inar", "hinar", "nbinar"),
                       nu_measles = NA, alpha_measles = NA, placeholder_measles = NA, psi_measles = NA, AIC_measles = NA,
                       nu_mumps = NA, alpha_mumps = NA, placeholder_mumps = NA, psi_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inar$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inar[i, "nu_measles"] <- fit_temp$coefficients["tau"]
  tab_inar[i, "alpha_measles"] <- fit_temp$coefficients["kappa"]
  tab_inar[i, "psi_measles"] <- fit_temp$coefficients["psi"]
  tab_inar[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inar[i, "nu_mumps"] <- fit_temp$coefficients["tau"]
  tab_inar[i, "alpha_mumps"] <- fit_temp$coefficients["kappa"]
  tab_inar[i, "psi_mumps"] <- fit_temp$coefficients["psi"] # WRONG
  tab_inar[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inar$model <- toupper(tab_inar$model)

c("Poisson INAR", "Hermite INAR", "NegBin INAR")



# INARCH:
# initialize table:
tab_inarch <- data.frame(model = c("inarch", "hinarch", "nbinarch"),
                       nu_measles = NA, alpha_measles = NA, placeholder_measles = NA, psi_measles = NA, AIC_measles = NA,
                       nu_mumps = NA, alpha_mumps = NA, placeholder_mumps = NA, psi_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inarch$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inarch[i, "nu_measles"] <- fit_temp$coefficients["nu"]
  tab_inarch[i, "alpha_measles"] <- fit_temp$coefficients["alpha"]
  tab_inarch[i, "psi_measles"] <- fit_temp$coefficients["psi"]  # WRONG
  tab_inarch[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inarch[i, "nu_mumps"] <- fit_temp$coefficients["nu"]
  tab_inarch[i, "alpha_mumps"] <- fit_temp$coefficients["alpha"]
  tab_inarch[i, "psi_mumps"] <- fit_temp$coefficients["psi"]  # WRONG
  tab_inarch[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inarch$model <- c("Poisson INARCH", "Hermite INARCH", "NegBin INARCH")



# INARMA:
# initialize table:
tab_inarma <- data.frame(model = c("inarma", "hinarma", "nbinarma"),
                       nu_measles = NA, kappa_measles = NA, beta_measles = NA, psi_measles = NA, AIC_measles = NA,
                       nu_mumps = NA, kappa_mumps = NA, beta_mumps = NA, psi_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inarma$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inarma[i, "nu_measles"] <- fit_temp$coefficients["tau"]
  tab_inarma[i, "kappa_measles"] <- fit_temp$coefficients["kappa"]
  tab_inarma[i, "beta_measles"] <- fit_temp$coefficients["beta"]
  tab_inarma[i, "psi_measles"] <- fit_temp$coefficients["psi"]
  tab_inarma[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inarma[i, "nu_mumps"] <- fit_temp$coefficients["tau"]
  tab_inarma[i, "kappa_mumps"] <- fit_temp$coefficients["kappa"]
  tab_inarma[i, "beta_mumps"] <- 1 - fit_temp$coefficients["phi"]
  tab_inarma[i, "psi_mumps"] <- fit_temp$coefficients["psi"]
  tab_inarma[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inarma$model <- c("Poisson INARMA", "Hermite INARMA", "NegBin INARMA")


# INGARCH:
# initialize table:
tab_ingarch <- data.frame(model = c("ingarch", "hingarch", "nbingarch"),
                         nu_measles = NA, alpha_measles = NA, beta_measles = NA, psi_measles = NA, AIC_measles = NA,
                         nu_mumps = NA, alpha_mumps = NA, beta_mumps = NA, psi_mumps = NA, AIC_mumps = NA)
# fille table:
for(i in 1:3){
  mod <- tab_ingarch$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_ingarch[i, "nu_measles"] <- fit_temp$coefficients["nu"]
  tab_ingarch[i, "alpha_measles"] <- fit_temp$coefficients["alpha"]
  tab_ingarch[i, "beta_measles"] <- fit_temp$coefficients["beta"]
  tab_ingarch[i, "psi_measles"] <- fit_temp$coefficients["psi"]
  tab_ingarch[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_ingarch[i, "nu_mumps"] <- fit_temp$coefficients["nu"]
  tab_ingarch[i, "alpha_mumps"] <- fit_temp$coefficients["alpha"]
  tab_ingarch[i, "beta_mumps"] <- fit_temp$coefficients["beta"]
  tab_ingarch[i, "psi_mumps"] <- fit_temp$coefficients["psi"]
  tab_ingarch[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_ingarch$model <- c("Poisson INGARCH", "Hermite INGARCH", "NegBin INGARCH")

# write out:
library(xtable)

for(tab in c("tab_inar", "tab_inarch", "tab_inarma", "tab_ingarch")){
  tab_temp <- get(tab)

  write(
    print(xtable(tab_temp), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("Tables/", tab, ".tex")
  )
}


#####################################
###### epidemiological scales:


# INAR:
# initialize table:
tab_inar_epi <- data.frame(model = c("inar", "hinar", "nbinar"),
                           immig_measles = NA, R_measles = NA, GT_measles = NA, CS_measles = NA, AIC_measles = NA,
                           immig_mumps = NA, R_mumps = NA,  GT_mumps = NA, CS_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inar_epi$model[i]

  fit_temp <- fits_measles[[mod]]
  tab_inar_epi[i, "R_measles"] <- fit_temp$coefficients["kappa"]
  tab_inar_epi[i, "immig_measles"] <- fit_temp$coefficients["tau"]
  # generation time is 1 by construction:
  tab_inar_epi[i, "GT_measles"] <- 1
  # cluster size is 1 by construction:
  tab_inar_epi[i, "CS_measles"] <- 1
  tab_inar_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inar_epi[i, "R_mumps"] <- fit_temp$coefficients["kappa"]
  tab_inar_epi[i, "immig_mumps"] <- fit_temp$coefficients["tau"]
  # generation time is 1 by construction:
  tab_inar_epi[i, "GT_mumps"] <- 1
  # cluster size is 1 by construction:
  tab_inar_epi[i, "CS_mumps"] <- 1
  tab_inar_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inar_epi$model <- c("Poisson INAR", "Hermite INAR", "NegBin INAR")



# INARCH:
# initialize table:
tab_inarch_epi <- data.frame(model = c("inarch", "hinarch", "nbinarch"),
                             immig_measles = NA, R_measles = NA, GT_measles = NA, CS_measles = NA, AIC_measles = NA,
                             immig_mumps = NA, R_mumps = NA, GT_mumps = NA, CS_mumps = NA, AIC_mumps = NA)
# fill table:
for(i in 1:3){
  mod <- tab_inarch_epi$model[i]
  family <- switch(mod, inarch = "Poisson", hinarch = "Hermite", nbinarch = "NegBin")

  fit_temp <- fits_measles[[mod]]
  tab_inarch_epi[i, "R_measles"] <- fit_temp$coefficients["alpha"]
  tab_inarch_epi[i, "immig_measles"] <- fit_temp$coefficients["nu"]
  # generation time is 1 by construction:
  tab_inarch_epi[i, "GT_measles"] <- 1
  # computing cluster size via custom function:
  tab_inarch_epi[i, "CS_measles"] <-
    get_cluster_size(psi = fit_temp$coefficients["psi"], family = family)
  tab_inarch_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_inarch_epi[i, "R_mumps"] <- fit_temp$coefficients["alpha"]
  tab_inarch_epi[i, "immig_mumps"] <- fit_temp$coefficients["nu"]
  # generation time is 1 by construction:
  tab_inarch_epi[i, "GT_mumps"] <- 1
  # computing cluster size via custon function:
  tab_inarch_epi[i, "CS_mumps"] <-
    get_cluster_size(psi = fit_temp$coefficients["psi"], family = family)
  tab_inarch_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_inarch_epi$model <- c("Poisson INARCH", "Hermite INARCH", "NegBin INARCH")



# INARMA:
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
  tab_ingarch_epi[i, "R_measles"] <- fit_temp$coefficients["alpha"]/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "immig_measles"] <- fit_temp$coefficients["nu"]/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "GT_measles"] <- 1/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "CS_measles"] <- get_cluster_size(psi = fit_temp$coefficients["psi"], family = family)
  tab_ingarch_epi[i, "AIC_measles"] <- fit_temp$AIC

  fit_temp <- fits_mumps[[mod]]
  tab_ingarch_epi[i, "R_mumps"] <- fit_temp$coefficients["alpha"]/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "immig_mumps"] <- fit_temp$coefficients["nu"]/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "GT_mumps"] <- 1/(1 - fit_temp$coefficients["beta"])
  tab_ingarch_epi[i, "CS_mumps"] <- get_cluster_size(psi = fit_temp$coefficients["psi"], family = family)
  tab_ingarch_epi[i, "AIC_mumps"] <- fit_temp$AIC
}
tab_ingarch_epi$model <- c("Poisson INGARCH", "Hermite INGARCH", "NegBin INGARCH")

# write out:
library(xtable)

for(tab in c("tab_inar_epi", "tab_inarch_epi", "tab_inarma_epi", "tab_ingarch_epi")){
  tab_temp <- get(tab)

  write(
    print(xtable(tab_temp), only.contents = TRUE,
          include.rownames = FALSE, include.colnames = FALSE,
          hline.after = NULL),
    file = paste0("Tables/", tab, ".tex")
  )
}


#####################################
###### diagnostic plot

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
pdf("Plots/residuals_inarma.pdf", width = 8, height = 3.5)

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
pdf("Plots/residuals_ingarch.pdf", width = 8, height = 3.5)

# define colours
cols_inarch <- brewer.pal(6, "Paired")
names(cols_inarch) <- c("inarch", "inarch",
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
