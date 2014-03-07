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

  ## summarizing function
  summarize <-
    function(mu){
      mean <- mu
      return(data.frame(mean = mean,
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(mu, n){
      return(rep(mu, n))
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
           summarize = summarize,
           sim = sim,
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

  ## summarizing function
  summarize <-
    function(mu, k){
      gamma_shape <- k
      gamma_rate  <- k / mu
      mean <- mu
      sdev <- sqrt(gamma_shape / (gamma_rate ^ 2))
      cnfi <- qgamma(c(0.025, 0.975), shape = gamma_shape, rate = gamma_rate)
      return(data.frame(mean = mean, sd = sdev,
                        "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(mu, k, n){
      return(rgamma(n, shape = k, rate = k / mu))
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
           summarize = summarize,
           sim = sim,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Poisson-Log Normal ------------------------------------------------------
poislognorm <-
function(x, q){
  ## summarizing function
  summarize <- lognormal()$summarize

  ## simulation function
  sim <- lognormal()$sim

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
           summarize = summarize,
           sim = sim,
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
    function(mu, shape, x, q){
      lik <-
        (exp(shape) / factorial(x)) *
        ((mu * q) ^ x) *
        ((shape * (shape + 2 * mu * q)) ^ (1 / 4)) *
        ((shape / (shape + 2 * mu * q)) ^ (x / 2)) *
        besselK(x = sqrt(shape * (shape + 2 * mu * q)), nu = x - 1 / 2) *
        sqrt(2 / pi)
      return(-sum(log(lik)))
    }

  ## summarizing function
  summarize <- invgauss()$summarize

  ## simulation function
  sim <- invgauss()$sim

  ## starting values
  start <- list(mu = 0.005, shape = 1)

  ## number of parameters
  npar <- 2

  ## return family
  return(
    structure(
      list(family = "poisinvgauss",
           minloglik = minloglik,
           summarize = summarize,
           sim = sim,
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

  ## summarizing function
  summarize <-
    function(eta, omega, lambda){
      Theta <- gigChangePars(4, 1, c(lambda, omega, eta))
      a <- Theta[3]
      b <- Theta[2]
      p <- Theta[1]
      mean <- sqrt(b) * besselK(sqrt(a * b), p + 1) /
             (sqrt(a) * besselK(sqrt(a * b), p))
      sdev <- sqrt((b / a) *
                   (besselK(sqrt(a * b), p + 2) / besselK(sqrt(a * b), p) -
                   (besselK(sqrt(a * b), p + 1) / besselK(sqrt(a * b), p)) ^ 2))
      cnfi <- qgig(c(.025, .975), Theta)
      return(data.frame(mean = mean, sd = sdev,
                        "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(eta, omega, lambda, n){
      return(rgig(n, gigChangePars(4, 1, c(lambda, omega, eta))))
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
           summarize = summarize,
           sim = sim,
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

  ## summarizing function
  summarize <-
    function(shape, rate){
      mean <- shape / rate
      sdev <- sqrt(shape) / rate
      cnfi <- qgamma(c(0.025, 0.975), shape, rate)
      return(data.frame(mean = mean, sd = sdev,
                        "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(shape, rate, n){
      return(rgamma(n, shape, rate))
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
           summarize = summarize,
           sim = sim,
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

  ## summarizing function
  summarize <-
    function(shape, scale){
      mean <- scale * base::gamma(1 + 1/shape)
      sdev <- scale * sqrt(base::gamma(1 + 2/shape) -
                           base::gamma(1 + 1/shape)^2)
      cnfi <- qweibull(c(0.025, 0.975), shape, scale)
      return(data.frame(mean = mean, sd = sdev,
                        "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(shape, scale, n){
      return(rweibull(n, shape, scale))
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
           summarize = summarize,
           sim = sim,
           start = start,
           npar = npar),
      class = "family")
  )
}

## Log-Normal --------------------------------------------------------------
lognormal <-
function(x = 1, d){
  ## minus log likelihood function
  minloglik <-
    function(mu_log, sd_log, x, d){
      d1 <- plnorm(x[d == 1], meanlog = mu_log, sdlog = sd_log)  # cens
      d0 <- dlnorm(x[d == 0], meanlog = mu_log, sdlog = sd_log)  # obs
      return(sum(-log(d1)) + sum(-log(d0)))
    }

  ## summarizing function
  summarize <-
    function(mu_log, sd_log){
      mean <- exp(mu_log + sd_log^2 / 2)
      sdev <- sqrt((exp(sd_log^2) - 1) * exp(2 * mu_log + sd_log^2))
      cnfi <- qlnorm(c(0.025, 0.975), mu_log, sd_log)
      return(data.frame(mean = mean, sd = sdev,
                        "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(mu_log, sd_log, n){
      return(rlnorm(n, mu_log, sd_log))
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
           summarize = summarize,
           sim = sim,
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

  ## summarizing function
  summarize <-
    function(mu, shape){
      mean <- mu
      sdev <- sqrt(mu^3 / shape)
      cnfi <- qinvGauss(c(0.025, 0.975), mu, shape)
      return(data.frame(mean = mean, sd = sdev,
                        "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                        check.names = FALSE, row.names = ""))
    }

  ## simulation function
  sim <-
    function(mu, shape, n){
      return(rinvGauss(n, mu, shape))
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
           summarize = summarize,
           sim = sim,
           start = start,
           npar = npar),
      class = "family")
  )
}