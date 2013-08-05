## FAMILIES ================================================================

## Poisson -----------------------------------------------------------------
poisson <-
function(x, q){
  ## minus log likelihood function (concentration)
  minloglik <-
    function(mu, x, q){
      lik <- (mu * q) ^ x * exp(-mu * q) / factorial(x)
      return(-sum(log(lik)))
    }

  ## minus log likelihood function (presence/absence)
  minloglik_bernoulli <-
    function(mu, x, q){
      lik <- (1 - exp(-mu * q)) ^ x * (exp(-mu * q)) ^ (1 - x)
      return(-sum(log(lik)))
    }

  ## starting values
  start <- list(mu = 1)

  ## number of parameters
  npar <- 1

  ## Bayesian implementation
  bayes <-
    list(prior = "mu ~ dgamma(0.01, 0.01)")

  ## return family
  return(
    structure(
      list(family = "poisson",
           minloglik = minloglik,
           minloglik_bernoulli = minloglik_bernoulli,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Negative Binomial -------------------------------------------------------
negbin <-
function(x, q){
  ## minus log likelihood function
  minloglik <-
    function(mu, k, x, q){
      lik <-
        (base::gamma(x + k) / (base::gamma(k) * factorial(x))) *
        ((mu * q / (k + mu * q)) ^ x) *
        (((k + mu * q) / k) ^ (-k))
      return(-sum(log(lik)))
    }

  ## starting values
  start <- list(mu = 0.005, k = 2)

  ## number of parameters
  npar <- 2

  ## Bayesian implementation
  bayes <-
    list(prior = c("mu ~ dgamma(shape, rate)",
                   "shape ~ dgamma(0.01, 0.01)",
                   "rate ~ dgamma(0.01, 0.01)"))

  ## return family
  return(
    structure(
      list(family = "negbin",
           minloglik = minloglik,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Poisson-Log Normal ------------------------------------------------------
poislognorm <-
function(x, q){
  ## integrand
  f <-
    function(mu, .x, .q, .mu_log, .sd_log)
      ((mu * .q) ^ .x * exp(-mu * .q) / factorial(.x)) *
      (1 / (mu * .sd_log * sqrt(2 * pi))) *
      (exp(-(log(mu) - .mu_log)^2 / (2 * .sd_log ^ 2)))

  ## minus log likelihood function
  minloglik <-
    function(mu_log, sd_log, x, q){
      lik <- numeric(length(x))
      for (i in seq_along(x)){
        lik[i] <-
          integrate(f, lower = 0, upper = Inf, .x = x[i], .q = q[i],
                    .mu_log = mu_log, .sd_log = sd_log,
                    stop.on.error = FALSE)$value
      }
      return(-sum(log(lik)))
    }

  ## starting values
  start <- list(mu_log = 0, sd_log = 1)

  ## number of parameters
  npar <- 2

  ## Bayesian implementation
  bayes <-
    list(prior = c("mu ~ dlnorm(mu_log, tau_log)",
                   "mu_log ~ dnorm(0, 0.0001)",
                   "tau_log ~ dgamma(0.01, 0.01)"))

  ## return family
  return(
    structure(
      list(family = "poislognorm",
           minloglik = minloglik,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Poisson-Inverse Gaussian ------------------------------------------------
poisinvgauss <-
function(x, q){
  ## minus log likelihood function
  minloglik <-
    function(mu, phi, x, q){
      lik <-
        (exp(phi) / factorial(x)) *
        ((mu * q) ^ x) *
        ((phi * (phi + 2 * mu * q)) ^ (1 / 4)) *
        ((phi / (phi + 2 * mu * q)) ^ (x / 2)) *
        besselK(x = sqrt(phi * (phi + 2 * mu * q)), nu = x - 1 / 2) *
        sqrt(2 / pi)
      return(-sum(log(lik)))
    }

  ## starting values
  start <- list(mu = 0.005, phi = 1)

  ## number of parameters
  npar <- 2

  ## return family
  return(
    structure(
      list(family = "poisinvgauss",
           minloglik = minloglik,
           start = start,
           npar = npar),
      class = "family")
  )
}

## Poisson-Generalized Inverse Gaussian ------------------------------------
poisgeninvgauss <-
function(x, q){
  ## minus log likelihood function
  minloglik <-
    function(eta, omega, lambda, x, q){
      lik <-
        (((eta * q) ^ x) / factorial(x)) *
        ((omega / (omega + 2 * eta * q)) ^ ((x + lambda) / 2)) *
        besselK(x = sqrt(omega * (omega + 2 * eta * q)), nu = lambda + x) /
        besselK(x = omega, nu = lambda)
      return(-sum(log(lik)))
    }

  ## number of parameters
  npar <- 3

  ## starting values
  start <- list(eta = 0.05, omega = 1, lambda = -1/2)

  ## return family
  return(
    structure(
      list(family = "poisgeninvgauss",
           minloglik = minloglik,
           start = start,
           npar = npar),
      class = "family")
  )
}

## Gamma -------------------------------------------------------------------
gamma <-
function(x, d){
  ## minus log likelihood function
  minloglik <-
    function(shape, rate, x, d){
      d1 <- pgamma(x[d == 1], shape = shape, rate = rate)  # cens
      d0 <- dgamma(x[d == 0], shape = shape, rate = rate)  # obs
      return(sum(-log(d1)) + sum(-log(d0)))
    }

  ## number of parameters
  npar <- 2

  ## starting values
  start <- list(shape = 1, rate = 1)

  ## return family
  return(
    structure(
      list(family = "gamma",
           minloglik = minloglik,
           start = start,
           npar = npar),
      class = "family")
  )
}

## Weibull -----------------------------------------------------------------
weibull <-
function(x, d){
  ## minus log likelihood function
  minloglik <-
    function(shape, scale, x, d){
      d1 <- pweibull(x[d == 1], shape = shape, scale = scale)  # cens
      d0 <- dweibull(x[d == 0], shape = shape, scale = scale)  # obs
      return(sum(-log(d1)) + sum(-log(d0)))
    }

  ## number of parameters
  npar <- 2

  ## starting values
  start <- list(shape = 1, scale = 1)

  ## return family
  return(
    structure(
      list(family = "weibull",
           minloglik = minloglik,
           start = start,
           npar = npar),
      class = "family")
  )
}

## Log-Normal --------------------------------------------------------------
lognormal <-
function(x, d){
  ## minus log likelihood function
  minloglik <-
    function(mu_log, sd_log, x, d){
      d1 <- plnorm(x[d == 1], meanlog = mu_log, sdlog = sd_log)  # cens
      d0 <- dlnorm(x[d == 0], meanlog = mu_log, sdlog = sd_log)  # obs
      return(sum(-log(d1)) + sum(-log(d0)))
    }

  ## number of parameters
  npar <- 2

  ## starting values
  start <- list(mu_log = mean(log(x)), sd_log = sd(log(x)))

  ## return family
  return(
    structure(
      list(family = "lognormal",
           minloglik = minloglik,
           start = start,
           npar = npar),
      class = "family")
  )
}

## Inverse Gaussian ---------------------------------------------------------
invgauss <-
function(x, d){
  ## minus log likelihood function
  minloglik <-
    function(mu, shape, x, d){
      y1 <- x[d == 1]  # cens
      d1 <- pnorm(((y1 / mu) - 1) * sqrt(mu * shape / y1)) +
            exp(2 * shape) *
            pnorm(-((y1 / mu) + 1) * sqrt(mu * shape / y1))
      y0 <- x[d == 0]  # obs
      d0 <- sqrt(mu * shape / (2 * pi * y0 ^ 3)) *
            exp(shape * (1 - 0.5 * ((y0 / mu) + (mu / y0))))
      return(sum(-log(d1)) + sum(-log(d0)))
    }

  ## number of parameters
  npar <- 2

  ## starting values
  start <- list(mu = 1, shape = 1)

  ## return family
  return(
    structure(
      list(family = "invgauss",
           minloglik = minloglik,
           start = start,
           npar = npar),
      class = "family")
  )
}