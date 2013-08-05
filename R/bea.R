## BAYESIAN EXPOSURE ASSESSMENT

## COUNT DATA --------------------------------------------------------------
bea_count <-
function(x, q = 1, data,
         family = c("poisson", "negbin", "poislognorm",
                    "poisinvgauss", "poisgeninvgauss"),
         nchains = 2, burnin = 1000, update = 5000, verbose = FALSE){

  ## check data
  if (missing(data) || is.null(data)){
    data <- parent.frame()

  } else if (!is.data.frame(data)){
    data <- data.frame(data)
  }

  ## match 'x' and 'q' arguments
  call <- match.call()
  call_x <- call[[match("x", names(call))]]
  call_q <- call[[match("q", names(call))]]
  x <- eval(call_x, data, enclos = parent.frame())
  q <- eval(call_q, data, enclos = parent.frame())

  ## check family
  family <- match.arg(family)
  family <- getFromNamespace(family, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N){",
      "x[i] ~ dpois(lambda[i])",
      "lambda[i] <- mu * q[i]",
      "}",
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(x = x, q = q, N = length(x))

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose)
    cat("JAGS progress:\n\n")

  JAGSout <-
    R2JAGS(model = model, data = data, inits = inits,
           nchains = nchains, burnin = burnin, update = update,
           nodes = "mu", verbose = verbose)

  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")

  DIC <- JAGSout$dic
  BGR <- c(gelman.diag(mcmc.list, autoburnin = FALSE)$psrf)

  bea_fit <-
    new("bea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        par = list(nchains = nchains, burnin = burnin,
                   update = update, inits = inits),
        model = model,
        mcmc = mcmc.list,
        diagnostics =
          list(DIC = DIC,
               BGR = data.frame(mean = BGR[1], upperCL = BGR[2])))

  ## return 'bea' object
  return(bea_fit)
}

## CONCENTRATION DATA ------------------------------------------------------
bea_conc <-
function(x, d, data,
         family = c("gamma", "lognormal", "weibull", "invgauss")){

  ## check data
  if (missing(data) || is.null(data)){
    data <- parent.frame()

  } else if (!is.data.frame(data)){
    data <- data.frame(data)
  }

  ## match 'x' and 'q' arguments
  call <- match.call()
  call_x <- call[[match("x", names(call))]]
  call_d <- call[[match("d", names(call))]]
  x <- eval(call_x, data, enclos = parent.frame())
  d <- eval(call_d, data, enclos = parent.frame())

  ## check family
  family <- match.arg(family)
  family <- getFromNamespace(family, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N){",
      "x[i] ~ dbern(pi[i])",
      "pi[i] <- 1 - exp(-lambda[i])",
      family()$bayes$model,
      "}",
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(x = x, q = q, N = length(x))

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose)
    cat("JAGS progress:\n\n")

  JAGSout <-
    R2JAGS(model = model, data = data, inits = inits,
           nchains = nchains, burnin = burnin, update = update,
           nodes = "mu", verbose = verbose)

  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")

  DIC <- JAGSout$dic
  BGR <- c(gelman.diag(mcmc.list, autoburnin = FALSE)$psrf)

  bea_fit <-
    new("bea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        par = list(nchains = nchains, burnin = burnin,
                   update = update, inits = inits),
        model = model,
        mcmc = mcmc.list,
        diagnostics =
          list(DIC = DIC,
               BGR = data.frame(mean = BGR[1], upperCL = BGR[2])))

  ## return 'bea' object
  return(bea_fit)
}

## PRESENCE/ABSENCE DATA ---------------------------------------------------
bea_presence <-
function(x, q = 1, replicates = rep(1, length(x)), data,
         family = c("poisson"),
         nchains = 2, burnin = 1000, update = 5000, verbose = FALSE){

  ## check data
  if (missing(data) || is.null(data)){
    data <- parent.frame()

  } else if (!is.data.frame(data)){
    data <- data.frame(data)
  }

  ## match 'x' and 'q' arguments
  call <- match.call()
  call_x <- call[[match("x", names(call))]]
  call_q <- call[[match("q", names(call))]]
  call_r <- call[[match("replicates", names(call))]]
  x <- eval(call_x, data, enclos = parent.frame())
  q <- eval(call_q, data, enclos = parent.frame())
  r <- eval(call_r, data, enclos = parent.frame())

  ## convert 'x' and 'q' if needed
  if (!identical(r, rep(1, length(x)))){
    reps <- rbind(x, r - x)
    x <- rep(rep(c(1, 0), length(x)), c(reps))
    q <- rep(q, r)
  }

  ## check 'family'
  family <- match.arg(family)
  family <- getFromNamespace(family, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N){",
      "x[i] ~ dbern(pi[i])",
      "pi[i] <- 1 - exp(-lambda[i])",
      "lambda[i] <- mu * q[i]",
      "}",
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(x = x, q = q, N = length(x))

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose)
    cat("JAGS progress:\n\n")

  JAGSout <-
    R2JAGS(model = model, data = data, inits = inits,
           nchains = nchains, burnin = burnin, update = update,
           nodes = "mu", verbose = verbose)

  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")

  DIC <- JAGSout$dic
  BGR <- c(gelman.diag(mcmc.list, autoburnin = FALSE)$psrf)

  bea_fit <-
    new("bea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        par = list(nchains = nchains, burnin = burnin,
                   update = update, inits = inits),
        model = model,
        mcmc = mcmc.list,
        diagnostics =
          list(DIC = DIC,
               BGR = data.frame(mean = BGR[1], upperCL = BGR[2])))

  ## return 'bea' object
  return(bea_fit)
}