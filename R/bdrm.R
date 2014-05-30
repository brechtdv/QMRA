## BAYESIAN DOSE-RESPONSE MODELLING

bdrm <-
function(x, n, dose, data,
         model = c("betapoisson", "bp",
                   "exponential", "exp",
                   "loglogistic", "ll",
                   "logprobit", "lp",
                   "extremevalue", "ev"),
         inits = NULL, nchains = 2, burnin = 10000, update = 10000,
         verbose = FALSE) {

  ## check data
  if (missing(data) || is.null(data)){
    data <- parent.frame()

  } else if (!is.data.frame(data)){
    data <- data.frame(data)
  }

  ## match data arguments
  call <- match.call()
  call_x <- call[[match("x", names(call))]]
  call_n <- call[[match("n", names(call))]]
  call_d <- call[[match("dose", names(call))]]
  x <- eval(call_x, data, enclos = parent.frame())
  n <- eval(call_n, data, enclos = parent.frame())
  d <- eval(call_d, data, enclos = parent.frame())

  ## check 'model'
  model <- match.arg(model)
  full <- c("betapoisson", "exponential",
            "loglogistic", "logprobit", "extremevalue")
  short <- c("bp", "exp", "ll", "lp", "ev")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## create Bayesian model
  model <-
    c("model {",
      "for (i in 1:N) {",
      "x[i] ~ dbin(p[i], n[i])",
      family()$bayes$model,
      "}",
      family()$bayes$prior,
      "}")

  ## define S3 class
  class(model) <- "JAGS_model"

  ## create data
  data <- list(N = length(x), x = x, n = n, dose = d)

  ## get results!
  if (verbose)
    cat("JAGS progress:\n\n")

  JAGSout <-
    R2JAGS(model = model, data = data, inits = inits,
           nchains = nchains, burnin = burnin, update = update,
           nodes = family()$bayes$nodes, verbose = verbose)

  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")

  DIC <- JAGSout$dic
  BGR <- c(gelman.diag(mcmc.list, autoburnin = FALSE)$psrf)

  bdrm_fit <-
    new("bdrm",
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

  ## return 'bdrm' object
  return(bdrm_fit)
}
