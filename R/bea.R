## BAYESIAN EXPOSURE ASSESSMENT

## COUNT DATA --------------------------------------------------------------
bea_count <-
function(x, q = 1, data,
         model = c("poisson", "p",
                   "negbin", "nb",
                   "poislognorm", "pln",
                   "poisweibull", "pw"),
         inits = NULL, nchains = 2, burnin = 5000, update = 5000,
         verbose = FALSE) {

  ## check data
  if (missing(data) || is.null(data)) {
    data <- parent.frame()

  } else if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  ## match 'x' and 'q' arguments
  call <- match.call()
  call_x <- call[[match("x", names(call))]]
  call_q <- call[[match("q", names(call))]]
  x <- eval(call_x, data, enclos = parent.frame())
  q <- eval(call_q, data, enclos = parent.frame())

  ## check 'model'
  model <- match.arg(model)
  full <- c("poisson", "negbin", "poislognorm", "poisweibull")
  short <- c("p", "nb", "pln", "pw")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N) {",
      "x[i] ~ dpois(lambda[i])",
      "lambda[i] <- mu * q[i]",
      "}",
      ifelse(is.null(family()$bayes$likelihood),
             "", paste("mu ~", family()$bayes$likelihood)),
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(N = length(x), x = x, q = q)

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
        data = as.data.frame(data[-1]),
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
         model = c("gamma", "g",
                   "lognorm", "ln",
                   "weibull", "w",
                   "invgauss", "ig"),
         inits = NULL, nchains = 2, burnin = 5000, update = 5000,
         verbose = FALSE) {

  ## check data
  if (missing(data) || is.null(data)) {
    data <- parent.frame()

  } else if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  ## match 'x' and 'q' arguments
  call <- match.call()
  call_x <- call[[match("x", names(call))]]
  call_d <- call[[match("d", names(call))]]
  x <- eval(call_x, data, enclos = parent.frame())
  d <- eval(call_d, data, enclos = parent.frame())

  ## redefine 'd' -> 0 == left-censored, 1 == observed
  d <- 1 - d

  ## define 'lod' -> 'lod' of observations must be < observation
  lod <- x
  lod[d == 1] <- 0

  ## redefine 'x' -> left-censored observations must be 'NA'
  x[d == 0] <- NA

  ## check 'model'
  model <- match.arg(model)
  full <- c("gamma", "lognorm", "weibull", "invgauss")
  short <- c("g", "ln", "w", "ig")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N) {",
      "d[i] ~ dinterval(x[i], lod[i])",
      paste("x[i] ~", family()$bayes$likelihood),
      "}",
      paste("mu ~", family()$bayes$likelihood),
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(N = length(x), x = x, d = d, lod = lod)

  ## generate inits
  ## left-censored observations must be initialized below 'lod'
  ## + inits must be > 0, for log-normal model
  if (is.null(inits)) {
    x_init <- rep(NA, length(x))
    x_init[d == 0] <- .Machine$double.eps
    inits <- list(x = x_init)
  }

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
        data = as.data.frame(data[-1]),
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
         model = c("poisson", "p"),
         inits = NULL, nchains = 2, burnin = 5000, update = 5000,
         verbose = FALSE) {

  ## check data
  if (missing(data) || is.null(data)) {
    data <- parent.frame()

  } else if (!is.data.frame(data)) {
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

  ## check 'model'
  model <- match.arg(model)
  full <- c("poisson")
  short <- c("p")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N) {",
      "x[i] ~ dbern(pi[i])",
      "pi[i] <- 1 - exp(-lambda[i])",
      "lambda[i] <- mu * q[i]",
      "}",
      ifelse(is.null(family()$bayes$likelihood),
             "", paste("mu ~", family()$bayes$likelihood)),
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(N = length(x), x = x, q = q)

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
        data = as.data.frame(data[-1]),
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
