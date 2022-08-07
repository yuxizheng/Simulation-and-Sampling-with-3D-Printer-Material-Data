#' Yuxi Zheng, s2265133, yuxizheng
#' Add your own function definitions at the top of this file.

#' Confidence Interval with Approximation
#'
#' Compute the confidence interval with method of approximation, given i.i.d. observations
#'
#' @param y vector of observations
#' @param alpha 1-alpha is the intended coverage probability
#' 
#' @return CI the confidence interval with 90% coverage
CI_approx <- function(y, alpha = 0.05) {
  
  n <- length(y)
  lambda_hat <- mean(y)
  
  # interval for given normal distribution
  theta_interval <-
    lambda_hat - sqrt(lambda_hat / n) * qnorm(c(1 - alpha / 2, alpha / 2))
  
  # make sure the interval is positive
  CI <- pmax(theta_interval, 0)
  
  CI
}

#' Confidence Interval with Non Approximation
#' 
#'  Compute the confidence interval with method of non-approximation, given i.i.d. observations
#' 
#' @param y vector of observations
#' @param alpha 1-alpha is the intended coverage probability
#' 
#' @return CI the confidence interval with 90% coverage
CI_non_approx <- function(y, alpha = 0.05) {
  
  n <- length(y)
  lambda_hat <- mean(y)
  a <- qnorm(alpha/2)
  b <- qnorm(1-alpha/2)
  
  # compute lower and upper
  lower <- b*b/(2*n) - b*sqrt(lambda_hat/n + b^2/(4*n*n)) + lambda_hat
  upper <- a*a/(2*n) - a*sqrt(lambda_hat/n + a^2/(4*n*n)) + lambda_hat
  theta_interval <- c(lower, upper)
  
  CI <- pmax(theta_interval, 0)
  
  CI
}

#' Log prior density
#' 
#' Compute the logarithm of the joint prior density
#' 
#' @param theta vector of theta parameter
#' @param params vector of gamma parameter
#' 
#' @return log_pd computed log prior density
log_prior_density <- function(theta, params) {
  # compute joint log density
  log_pd <- dnorm(theta[1], mean = 0, sd = sqrt(params[1]), log = TRUE) + 
  dnorm(theta[2], mean = 1, sd = sqrt(params[2]), log = TRUE) +
  dlogexp(theta[3], rate = params[3], log = TRUE) +
  dlogexp(theta[4], rate = params[4], log = TRUE)
  
  log_pd
}

#' Log Observation likelihood 
#' 
#' Compute the logarithm of the likelihood
#' 
#' @param theta vector of theta parameter
#' @param x vector of x
#' @param y vector of y
#' 
#' @return log_likelihood computed log likelihood 
log_like <- function(theta, x, y) {
  # compute log density for samples
  log_f <- dnorm(y, mean = theta[1] + theta[2] * x, sd = sqrt(exp(theta[3]) + exp(theta[4]) * x^2), log = TRUE)
  log_likelihood <- sum(log_f)
  
  log_likelihood 
}

#' Posterior density
#' 
#' Compute the log posterior density, ignoring un-evaluated part
#' 
#' @param theta vector of theta parameter
#' @param x vector of x
#' @param y vector of y
#' @param params vector of gamma parameter
#' 
#' @return log_post log_posterior density
log_posterior_density <- function(theta, x, y, params) {
  
  log_post <- log_prior_density(theta, params) + log_like(theta, x, y)
  
  log_post
}

#' Importance Sampling
#' 
#' Compute 90% Bayesian credible interval for each βj using importance sampling
#' 
#' @param N The number of samples to generate
#' @param mu The mean vector for the importance distribution
#' @param S The covariance matrix
#' 
#' @return df data frame of samples and log-weights
do_importance <- function(N, mu, S, x, y, params, ...) {
  
  sample <- rmvnorm(N, mean = mu, sigma = S)
  log_weights <- vector(length = N)
  
  # for each of the simulated theta, calculate log weights
  for (i in 1:N) {
    log_weights[i] <- log_posterior_density(theta = sample[i,], x = x, y = y, params = params) - 
      dmvnorm(sample[i,], mean = mu, sigma = S, log = TRUE)
  }
  
  # normalization
  log_weights <- log_weights - log_sum_exp(log_weights)
  
  # compute weights
  sample[, 3] <- exp(sample[, 3])
  sample[, 4] <- exp(sample[, 4])
  
  df <- data.frame(cbind(sample, log_weights))
  colnames(df) <- c("beta1", "beta2", "beta3", "beta4", "log_weights")
  
  df
}

#' Construct credible intervals using importance sampling
#' 
#' Compute credible intervals as dataframe 
#' 
#' @param x vector of samples
#' @param weights vector of weights
#' @param prob intended coverage probability
#' 
#' @param df data frame of confidence interval for samples
make_CI <- function(x, weights, prob) {
  
  alpha <- 1-prob
  theta_interval <- wquantile(x, probs = c(alpha/2, 1-alpha/2), weights = weights)
  
  df <- data.frame(t(theta_interval))
  colnames(df) <- c("Lower", "Upper")
  
  df
}

#' Construct confidence intervals for input data of 
#' 
#' Compute mean, sd, lwr, upr for the normal distribution model as list
#' 
#' @param x vector of CAD weights
#' @param theta vector of theta parameters
#' @param prob confidence level
#' 
#' @return model_pred data frame of prediction information
model_prediction <- function(x, theta, prob){
  
  n <- length(x)
  alpha <- (1-prob) / 2 
  mean <- vector()
  sd <- vector()
  lwr <- vector()
  upr <- vector()
  
  # compute upper and lower, mean and sd for the point
  for(i in 1:n){
    mean <- append(mean, theta[1] + theta[2] * x[i])
    sd <- append(sd, sqrt(exp(theta[3]) + exp(theta[4]) * x[i]^2))
    lwr <- append(lwr, mean[i] + qnorm(alpha) * sd[i])
    upr <- append(upr, mean[i] - qnorm(alpha) * sd[i])
  }
  model_pred <- list(mean, sd, lwr, upr)
  names(model_pred) <- c("mean", "sd", "lwr", "upr")
  
  model_pred
}

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}

#' Interval coverage estimation
#'
#' @param CI_method a function object taking arguments x and alpha,
#'   returning a 2-element vector.
#' @param N The number of simulation replications to use for the coverage estimate
#' @param alpha 1-alpha is the intended coverage probability
#' @param n The sample size
#' @param lambda The true lambda values for the Poisson(lambda) model
#'
#' @return A scalar between 0 and 1, estimating the coverage probability
#' for the `CI_method` intervals for the given parameters.

estimate_coverage <- function(CI_method,
                              N = 10000,
                              alpha = 0.1,
                              n = 2,
                              lambda = 3) {
  cover <- 0
  for (loop in seq_len(N)) {
    y <- rpois(n, lambda)
    ci <- CI_method(y, alpha)
    cover <- cover + ((ci[1] <= lambda) && (lambda <= ci[2]))
  }
  cover / N
}

