## DOSE RESPONSE MODEL FAMILIES ============================================

## Beta-Poisson ------------------------------------------------------------
betapoisson <-
function(x, n, dose){
  ## minus log likelihood function
  minloglik <-
    function(logitU, V, x, n, d) {
      ## alpha, beta > 0
      ## U ~ (0,1)
      U <- exp(logitU) / (1 + exp(logitU))
      alpha <- U * exp(V)
      beta <- (1 - U) * exp(V)
      p <- 1 - (1 + (d / beta)) ^ (-alpha)
      loglik <- x * log(p) + (n - x) * log(1 - p)
      return(-sum(loglik))
    }

  ## predict function
  predict <-
    function(coef, vcov, dose) {
      ## point estimate
      logitU <- coef[1]
      V      <- coef[2]
      U <- exp(logitU) / (1 + exp(logitU))
      alpha <- U * exp(V)
      beta <- (1 - U) * exp(V)
      p <- 1 - (1 + (dose / beta)) ^ (-alpha)

      ## standard error ~ delta method
      f <- paste("~ 1 - (1 + (", dose, " / (1 - (exp(x1) / (1 + exp(x1)))) * exp(x2))) ^
                 (-(exp(x1) / (1 + exp(x1))) * exp(x2))", sep = "")
      f <- sapply(as.list(f), as.formula)
      p_se <- delta_method(f, coef, vcov)

      return(c(p, p_se))
    }

  ## simulation function
  sim <-
    function(n, mu, Sigma, dose) {
      sim <- matrix(mvrnorm(n, mu, Sigma), nrow = n)
      logitU <- sim[, 1]
      V <- sim[, 2]
      U <- exp(logitU) / (1 + exp(logitU))
      alpha <- U * exp(V)
      beta <- (1 - U) * exp(V)

      dose <- matrix(dose, ncol = 1)
      p <- t(apply(dose, 1,
                   function(x) 1 - (1 + (x / beta)) ^ (-alpha)))
      if (is.matrix(p) && nrow(p) == 1) p <- c(p)
      return(p)
    }

  ## starting values
  start <- list(logitU = 1, V = 1)

  ## number of parameters
  npar <- 2

  ## Bayesian implementation
  bayes <-
    list(model = "p[i] <- 1 - (1 + (dose[i] / beta)) ^ (-alpha)",
         prior = c("alpha ~ dnorm(0, 1e-10)T(0,)",
                   "beta  ~ dnorm(0, 1e-10)T(0,)"),
         nodes = c("alpha", "beta"),
         sim =
           function(pars, dose) {
             p <- t(apply(dose, 1,
                          function(x) 1 - (1 + (x / pars[, 2])) ^ (-pars[, 1])))
             if (is.matrix(p) && nrow(p) == 1) p <- c(p)
             return(p)
           })

  ## return family
  return(
    structure(
      list(family = "betapoisson",
           minloglik = minloglik,
           predict = predict,
           sim = sim,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Exponential -------------------------------------------------------------
exponential <-
function(x, n, dose){
  ## minus log likelihood function
  minloglik <-
    function(theta, x, n, d){
      p <- 1 - exp(-d * theta)
      loglik <- x * log(p) + (n - x) * log(1 - p)
      return(-sum(loglik))
    }

  ## starting values
  start <- list(theta = 0.001)

  ## number of parameters
  npar <- 1

  ## Bayesian implementation
  bayes <-
    list(model = "p[i] <- 1 - exp(-dose[i] * theta)",
         prior = "theta ~ dnorm(0, 1e-5)T(0,)",
         nodes = "theta",
         sim =
           function(pars, dose) {
             p <- t(apply(dose, 1,
                          function(x) 1 - exp(-x * pars[, 1])))
             if (is.matrix(p) && nrow(p) == 1) p <- c(p)
             return(p)
           })

  ## return family
  return(
    structure(
      list(family = "exponential",
           minloglik = minloglik,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Log-Logistic ------------------------------------------------------------
loglogistic <-
function(x, n, dose){
  ## minus log likelihood function
  minloglik <-
    function(alpha, beta, x, n, d) {
      theta <- alpha + beta * log(d)
      p <- 1 / (1 + exp(-theta))
      loglik <- x * log(p) + (n - x) * log(1 - p)
      return(-sum(loglik))
    }

  ## predict function
  predict <-
    function(coef, vcov, dose){
      ## point estimate
      theta <- coef[1] + coef[2] * log(dose)
      p <- 1 / (1 + exp(-theta))

      ## standard error ~ delta method
      f <- paste0("~ 1 / (1 + exp(-(x1 + x2 * log(", dose, "))))")
      f <- sapply(as.list(f), as.formula)
      p_se <- delta_method(f, coef, vcov)

      return(c(p, p_se))
    }

  ## simulation function
  sim <-
    function(n, mu, Sigma, dose) {
      sim <- matrix(mvrnorm(n, mu, Sigma), nrow = n)
      alpha <- sim[, 1]
      beta  <- sim[, 2]

      dose <- matrix(dose, ncol = 1)
      p <- t(apply(dose, 1,
                   function(x) 1 / (1 + exp(-(alpha + beta * log(x))))))
      if (is.matrix(p) && nrow(p) == 1) p <- c(p)
      return(p)
    }

  ## starting values
  start <- list(alpha = -1, beta = 1)

  ## number of parameters
  npar <- 2

  ## Bayesian implementation
  bayes <-
    list(model = c("p[i] <- 1 / (1 + exp(-theta[i]))",
                   "theta[i] <- alpha + beta * log(dose[i])"),
         prior = c("alpha ~ dnorm(0, 1e-5)T(,0)",
                   "beta  ~ dnorm(0, 1e-5)T(0,)"),
         nodes = c("alpha", "beta"),
         sim =
           function(pars, dose) {
             p <- t(apply(dose, 1,
                          function(x) 1 / (1 + exp(-(pars[, 1] + pars[, 2] * log(x))))))
             if (is.matrix(p) && nrow(p) == 1) p <- c(p)
             return(p)
           })

  ## return family
  return(
    structure(
      list(family = "loglogistic",
           minloglik = minloglik,
           predict = predict,
           sim = sim,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Log-Probit --------------------------------------------------------------
logprobit <-
function(x, n, dose){
  ## minus log likelihood function
  minloglik <-
    function(alpha, beta, x, n, d) {
      theta <- alpha + beta * log(d)
      p <- pnorm(theta)
      loglik <- x * log(p) + (n - x) * log(1 - p)
      return(-sum(loglik))
    }

  ## predict function
  predict <-
    function(coef, vcov, dose) {
      theta <- coef[1] + coef[2] * log(dose)
      p <- pnorm(theta)

      ## standard error ~ delta method
      f <- paste0("~ pnorm(x1 + x2 * log(", dose, "))")
      f <- sapply(as.list(f), as.formula)
      p_se <- delta_method(f, coef, vcov)

      return(c(p, p_se))
    }

  ## simulation function
  sim <-
    function(n, mu, Sigma, dose) {
      sim <- matrix(mvrnorm(n, mu, Sigma), nrow = n)
      alpha <- sim[, 1]
      beta  <- sim[, 2]

      dose <- matrix(dose, ncol = 1)
      p <- t(apply(dose, 1, 
                   function(x) pnorm(alpha + beta * log(x))))
      if (is.matrix(p) && nrow(p) == 1) p <- c(p)
      return(p)
    }

  ## starting values
  start <- list(alpha = -2, beta = 0.1)

  ## number of parameters
  npar <- 2

  ## Bayesian implementation
  bayes <-
    list(model = c("p[i] <- phi(theta[i])",
                   "theta[i] <- alpha + beta * log(dose[i])"),
         prior = c("alpha ~ dnorm(0, 1e-5)T(,0)",
                   "beta  ~ dnorm(0, 1e-5)T(0,)"),
         nodes = c("alpha", "beta"),
         sim =
           function(pars, dose) {
             p <- t(apply(dose, 1,
                          function(x) pnorm(pars[, 1] + pars[, 2] * log(x))))
             if (is.matrix(p) && nrow(p) == 1) p <- c(p)
             return(p)
           })

  ## return family
  return(
    structure(
      list(family = "logprobit",
           minloglik = minloglik,
           predict = predict,
           sim = sim,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}

## Extreme Value -----------------------------------------------------------
extremevalue <-
function(x, n, dose){
  ## minus log likelihood function
  minloglik <-
    function(alpha, beta, x, n, d){
      theta <- alpha + beta * log(d)
      p <- 1 - exp(-exp(theta))
      loglik <- x * log(p) + (n - x) * log(1 - p)
      return(-sum(loglik))
    }

  ## predict function
  predict <-
    function(coef, vcov, dose){
      ## point estimate
      theta <- coef[1] + coef[2] * log(dose)
      p <- 1 - exp(-exp(theta))

      ## standard error ~ delta method
      f <- paste0("~ 1 - exp(-exp(x1 + x2 * log(", dose, ")))")
      f <- sapply(as.list(f), as.formula)
      p_se <- delta_method(f, coef, vcov)

      return(c(p, p_se))
    }

  ## simulation function
  sim <-
    function(n, mu, Sigma, dose) {
      sim <- matrix(mvrnorm(n, mu, Sigma), nrow = n)
      alpha <- sim[, 1]
      beta  <- sim[, 2]

      dose <- matrix(dose, ncol = 1)
      p <- t(apply(dose, 1,
                   function(x) 1 - exp(-exp(alpha + beta * log(x)))))
      if (is.matrix(p) && nrow(p) == 1) p <- c(p)
      return(p)
    }

  ## starting values
  start <- list(alpha = -2, beta = 0.1)

  ## number of parameters
  npar <- 2

  ## Bayesian implementation
  bayes <-
    list(model = c("p[i] <- 1 - exp(-exp(theta[i]))",
                   "theta[i] <- alpha + beta * log(dose[i])"),
         prior = c("alpha ~ dnorm(0, 1e-5)T(,0)",
                   "beta  ~ dnorm(0, 1e-5)T(0,)"),
         nodes = c("alpha", "beta"),
         sim =
           function(pars, dose) {
             p <- t(apply(dose, 1,
                          function(x) 1 - exp(-exp(pars[, 1] + pars[, 2] * log(x)))))
             if (is.matrix(p) && nrow(p) == 1) p <- c(p)
             return(p)
           })

  ## return family
  return(
    structure(
      list(family = "extremevalue",
           minloglik = minloglik,
           predict = predict,
           sim = sim,
           start = start,
           npar = npar,
           bayes = bayes),
      class = "family")
  )
}
