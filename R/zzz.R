##= Define S4 classes =====================================================
setClass("ea",
  representation(
    mle_call = "language",
    data = "data.frame",
    family = "function",
    estimates = "list",
    AIC = "numeric"),
  contains = "mle")

setOldClass("JAGS_model") # virtual S3 class
setOldClass("mcmc.list") # virtual S3 class

setClass("bea",
  representation(
    call = "language",
    data = "data.frame",
    family = "function",
    par = "list",
    model = "JAGS_model",
    mcmc = "list",
    diagnostics = "list"))

##= Define S4 methods =====================================================
setMethod("initialize", "ea",
  function(.Object, call, data, family, estimates, AIC, mle, ...) {
    .Object <- callNextMethod()

    .Object@call <- call
    .Object@data <- data
    .Object@family <- family
    .Object@estimates <- estimates
    .Object@AIC <- AIC

    .Object@mle_call  <- mle@call
    .Object@coef      <- mle@coef
    .Object@fullcoef  <- mle@fullcoef
    .Object@vcov      <- mle@vcov
    .Object@min       <- mle@min
    .Object@details   <- mle@details
    .Object@minuslogl <- mle@minuslogl
    .Object@nobs      <- mle@nobs
    .Object@method    <- mle@method

    return(.Object)
  }
)

setMethod("show", "ea",
  function(object)
    print(object)
)

setMethod("print", "ea",
  function(x, dig = 3, ...){
    from <-
      c("count", "concentration", "presence/absence")[
        which(x@call[[1]] == c("ea_count", "ea_conc", "ea_presence"))]
    cat("Exposure assessment from", from, "data\n\n")
    cat("Call:\n")
    print(x@call)
    cat("\nCoefficients:\n")
    print(x@coef)
    cat("\nAIC:", x@AIC)
    cat("\n\n")
  }
)

setMethod("show", "bea",
  function(object)
    print(object)
)

setMethod("print", "bea",
  function(x, dig = 3, ...){
    from <-
      c("count", "concentration", "presence/absence")[
        which(x@call[[1]] == c("bea_count", "bea_conc", "bea_presence"))]
    cat("Bayesian exposure assessment from", from, "data\n\n")

    cat("Call:\n")
    print(x@call)

    cat("\nEstimate:\n")
    est <-
      data.frame(mean(unlist(x@mcmc), na.rm = TRUE),
                 quantile(unlist(x@mcmc), .025, na.rm = TRUE),
                 quantile(unlist(x@mcmc), .975, na.rm = TRUE))
    colnames(est) <- c("mean", "2.5%", "97.5%")
    rownames(est) <- "mu"
    print(est, digits = dig, ...)

    cat("\nDeviance Information Criterion:\n")
    print(x@diagnostics$DIC)

    cat("\nBGR statistic: ",
        formatC(x@diagnostics$BGR[, 1], digits = dig, format = "f"),
        " (upperCL ",
        formatC(x@diagnostics$BGR[, 2], digits = dig, format = "f"),
        ")", sep = "")
    cat("\nBGR values significantly higher than one indicate lack of fit.")
    cat("\n\n")
  }
)