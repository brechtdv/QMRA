##= Define S4 classes =====================================================
setClass("ea",
  representation(
    mle_call = "language",
    data = "data.frame",
    family = "function",
    estimates = "list",
    AIC = "numeric"),
  contains = "mle")

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