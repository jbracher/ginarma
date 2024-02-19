# Evaluate log-likelihood for an INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, S1: parameters of the INGARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_ingarch <- function(vect, nu, alpha, beta, lambda1, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dpois(vect, lambda, log = TRUE))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Fitting an INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_ingarch <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_ingarch(vect, nu = exp(pars["log_nu"]),
                  alpha = exp(pars["log_alpha"]),
                  beta = exp(pars["log_beta"]),
                  lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, log_beta = -1, log_lambda1 = 0.5)
  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                    alpha = exp(opt$par["log_alpha"]),
                    beta = exp(opt$par["log_beta"]),
                    lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "beta", "lambda1")
  AIC <- 2*(-llik + length(coefficients))

  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1] +
      coefficients["beta"]*fitted_values[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/sqrt(fitted_values)

  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC, opt = opt,
              fitted_values = fitted_values,
              pearson_residuals = pearson_residuals))
}


fit_inarch <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_ingarch(vect, nu = exp(pars["log_nu"]),
                  alpha = exp(pars["log_alpha"]),
                  beta = 0,
                  lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, log_lambda1 = 0.5)
  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                                alpha = exp(opt$par["log_alpha"]),
                                lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "lambda1")
  AIC <- 2*(-llik + length(coefficients))

  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/sqrt(fitted_values)

  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC, opt = opt,
              fitted_values = fitted_values,
              pearson_residuals = pearson_residuals))
}


# Evaluate log-likelihood for a negative binomial INGARCH(1, 1) model, version 1: constant size parameter
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, psi, S1: parameters of the INGARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_nbingarch <- function(vect, nu, alpha, beta, psi, lambda1, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dnbinom(vect, mu = lambda, size = 1/psi, log = TRUE))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Fitting a negative binomial INGARCH(1, 1) model, version 1: constant size parameter
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_nbingarch <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_nbingarch(vect, nu = exp(pars["log_nu"]),
                    alpha = exp(pars["log_alpha"]),
                    beta = exp(pars["log_beta"]),
                    psi = exp(pars["log_psi"]),
                    lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, log_beta = -1, log_psi = -0.5, log_lambda1 = 0.5)

  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                                alpha = exp(opt$par["log_alpha"]),
                                beta = exp(opt$par["log_beta"]),
                                psi = exp(opt$par["log_psi"]),
                                lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "beta", "psi", "lambda1")
  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1] +
      coefficients["beta"]*fitted_values[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/
    sqrt(fitted_values + coefficients["psi"]*fitted_values^2)

  AIC <- 2*(-llik + length(coefficients))
  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC,
              fitted_values = fitted_values, observed = vect,
              pearson_residuals = pearson_residuals, opt = opt))
}



# Evaluate log-likelihood for a negative binomial INGARCH(1, 1) model, version 2: size proportional to mean
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, psi, S1: parameters of the INGARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_nbingarch2 <- function(vect, nu, alpha, beta, psi, lambda1, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dnbinom(vect, mu = lambda, size = lambda/psi, log = TRUE))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Fitting a negative binomial INGARCH(1, 1) model, version 2: size proportional to mean
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_nbingarch2 <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_nbingarch2(vect, nu = exp(pars["log_nu"]),
                     alpha = exp(pars["log_alpha"]),
                     beta = exp(pars["log_beta"]),
                     psi = exp(pars["log_psi"]),
                     lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, log_beta = -1, log_psi = -0.5, log_lambda1 = 0.5)
  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                                alpha = exp(opt$par["log_alpha"]),
                                beta = exp(opt$par["log_beta"]),
                                psi = exp(opt$par["log_psi"]),
                                lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "beta", "psi", "lambda1")
  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1] +
      coefficients["beta"]*fitted_values[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/
    sqrt(fitted_values + coefficients["psi"]*fitted_values)

  AIC <- 2*(-llik + length(coefficients))
  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC,
              fitted_values = fitted_values, observed = vect,
              pearson_residuals = pearson_residuals, opt = opt))
}



# Fitting a negative binomial INARCH(1) model, version 2: size proportional to mean
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_nbinarch2 <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_nbingarch2(vect, nu = exp(pars["log_nu"]),
                     alpha = exp(pars["log_alpha"]),
                     beta = 0,
                     psi = exp(pars["log_psi"]),
                     lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, log_psi = -0.5, log_lambda1 = 0.5)
  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                                alpha = exp(opt$par["log_alpha"]),
                                psi = exp(opt$par["log_psi"]),
                                lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "psi", "lambda1")
  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/
    sqrt(fitted_values + coefficients["psi"]*fitted_values)

  AIC <- 2*(-llik + length(coefficients))
  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC,
              fitted_values = fitted_values, observed = vect,
              pearson_residuals = pearson_residuals, opt = opt))
}


# Evaluate log-likelihood for a Hermite INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, psi, S1: parameters of the INGARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_hingarch <- function(vect, nu, alpha, beta, psi, lambda1, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(inarma:::dherm(vect, mu = lambda, psi = psi, log = TRUE))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Fitting a Hermite INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_hingarch <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_hingarch(vect, nu = exp(pars["log_nu"]),
                   alpha = exp(pars["log_alpha"]),
                   beta = exp(pars["log_beta"]),
                   psi = exp(pars["logit_psi"])/(1 + exp(pars["logit_psi"])),
                   lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, log_beta = -1, logit_psi = -0.5, log_lambda1 = 0.5)
  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                                alpha = exp(opt$par["log_alpha"]),
                                beta = exp(opt$par["log_beta"]),
                                psi = exp(opt$par["logit_psi"]) / (1 + exp(opt$par["logit_psi"])),
                                lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "beta", "psi", "lambda1")
  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1] +
      coefficients["beta"]*fitted_values[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/
    sqrt((1 + coefficients["psi"])*fitted_values)

  AIC <- 2*(-llik + length(coefficients))
  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC,
              fitted_values = fitted_values, observed = vect,
              pearson_residuals = pearson_residuals, opt = opt))
}


fit_hinarch <- function(vect, start = NULL, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_hingarch(vect, nu = exp(pars["log_nu"]),
                   alpha = exp(pars["log_alpha"]),
                   beta = 0,
                   psi = exp(pars["logit_psi"])/(1 + exp(pars["logit_psi"])),
                   lambda1 = exp(pars["log_lambda1"]))
  }
  if(is.null(start)) start <- c(log_nu = 2, log_alpha = -1, logit_psi = -0.5, log_lambda1 = 0.5)
  opt <- optim(start, nllik,...)
  # structure results
  llik <- -opt$value
  coefficients_raw <- opt$par
  coefficients <- c(nu = exp(opt$par["log_nu"]),
                                alpha = exp(opt$par["log_alpha"]),
                                psi = exp(opt$par["logit_psi"])/(1 + exp(opt$par["logit_psi"])),
                                lambda1 = exp(opt$par["log_lambda1"]))
  names(coefficients) <- c("nu", "alpha", "psi", "lambda1")
  # compute fitted values:
  fitted_values <- numeric(length(vect))
  fitted_values[1] <- coefficients["lambda1"]
  for(i in 2:length(fitted_values)){
    fitted_values[i] <- coefficients["nu"] +
      coefficients["alpha"]*vect[i - 1]
  }
  pearson_residuals <- (vect - fitted_values)/
    sqrt((1 + coefficients["psi"])*fitted_values)

  AIC <- 2*(-llik + length(coefficients))
  return(list(coefficients = coefficients,
              coefficients_raw = coefficients_raw,
              llik = llik, dim = length(coefficients), AIC = AIC,
              fitted_values = fitted_values, observed = vect,
              pearson_residuals = pearson_residuals, opt = opt))
}


my_lines <- function(x, y1, y2, col1 = 1, col2 = 2, lwd = 2){
  inds1larger <- which(abs(y1) > abs(y2))
  inds2larger <- which(abs(y2) > abs(y1))
  lines(x[inds1larger], y1[inds1larger], type = "h", col = col1, lwd = lwd)
  lines(x[inds2larger], y2[inds2larger], type = "h", col = col2, lwd = lwd)
  lines(x[inds2larger], y1[inds2larger], type = "h", col = col1, lwd = lwd)
  lines(x[inds1larger], y2[inds1larger], type = "h", col = col2, lwd = lwd)
}

# compute the clustersize for the summary table from the model parameter psi
# and the distribuiton family
get_cluster_size <- function(psi, family = c("Poisson", "Hermite", "NegBin")){
  # computing cluster size is tedious:
  if(family == "Poisson"){
    clustersize <- 1
  }
  if(family == "Hermite"){
    pi <- psi / (2 - psi)
    clustersize <- 1 + pi
  }
  if(family == "NegBin"){
    pi <- 1/(1 + psi)
    clustersize <- (1 - pi)/(-pi*log(pi))
  }
  return(clustersize)
}
