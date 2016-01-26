##= Define S4 classes =====================================================

## first define some virtual S3 classes
setOldClass("JAGS_model")
setOldClass("mcmc.list")
setOldClass("gof") 

## 'ea': exposure assessment
setClass("ea",
  representation(
    mle_call = "language",
    data = "data.frame",
    family = "function",
    gof = "gof",
    AIC = "numeric",
    mle = "mle"),
  contains = "mle")

## 'bea': Bayesian exposure assessment
setClass("bea",
  representation(
    call = "language",
    data = "data.frame",
    family = "function",
    par = "list",
    model = "JAGS_model",
    mcmc = "list",
    diagnostics = "list"))

## 'drm': dose-response modelling
setClass("drm",
  representation(
    mle_call = "language",
    data = "data.frame",
    family = "function",
    AIC = "numeric",
    gof = "gof",
    mle = "mle"),
  contains = "mle")

## 'bdrm': Bayesian dose-response modelling
setClass("bdrm",
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
  function(.Object, call, data, family, gof, AIC, mle, ...) {
    .Object <- callNextMethod()

    .Object@call <- call
    .Object@data <- data
    .Object@family <- family
    .Object@gof <- gof
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

setMethod("initialize", "drm",
  function(.Object, call, data, family, AIC, gof, mle, ...) {
    .Object <- callNextMethod()

    .Object@call <- call
    .Object@data <- data
    .Object@family <- family
    .Object@AIC <- AIC
    .Object@gof <- gof

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
    getMethod("print", "ea")(object)
)

setMethod("print", "ea",
  function(x, dig = 3, ...) {
    ## data type
    from <-
      c("count", "concentration", "presence/absence")[
        which(x@call[[1]] == c("ea_count", "ea_conc", "ea_presence"))]
    cat("Exposure assessment from", from, "data\n\n")

    ## function call
    cat("Call:\n")
    print(x@call)

    ## model family & estimates
    family <-
      switch(x@family()$family,
             gamma = "Gamma",
             lognormal = "Log-normal",
             weibull = "Weibull",
             invgaus = "Inverse Gaussian",
             poisson = "Poisson",
             negbin = "Negative Binomial",
             poislognorm = "Poisson-Log-Normal",
             poisinvgauss = "Poisson-Inverse Gaussian",
             poisgeninvgauss = "Poisson-Generalised Inverse Gaussian")
    cat("\n", family, " model coefficients:\n", sep = "")
    print(x@coef)

    ## estimates
    cat("\nEstimated concentration:\n")
    print(summarize(x))

    ## AIC
    cat("\nAIC:", x@AIC)
    cat("\n\n")
  }
)

setMethod("show", "bea",
  function(object)
    getMethod("print", "bea")(object)
)

setMethod("print", "bea",
  function(x, dig = 3, ...) {
    from <-
      c("count", "concentration", "presence/absence")[
        which(x@call[[1]] == c("bea_count", "bea_conc", "bea_presence"))]
    cat("Bayesian exposure assessment from", from, "data\n\n")

    cat("Call:\n")
    print(x@call)

    cat("\nEstimated concentration:\n")
    est <- summarize(x)
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

setMethod("show", "drm",
  function(object)
    getMethod("print", "drm")(object)
)

setMethod("print", "drm",
  function(x, dig = 3, ...) {
    ## title
    cat("Dose-response model\n\n")

    ## function call
    cat("Call:\n")
    print(x@call)

    ## model family & estimates
    family <-
      switch(x@family()$family,
             betapoisson = "Beta-Poisson",
             exponential = "Exponential",
             loglogistic = "Log-Logistic",
             logprobit = "Log-Probit",
             extremevalue = "Extreme value")
    cat("\n", family, " model coefficients:\n", sep = "")
    print(coef(summary(x)))

    ## AIC
    cat("\nAIC:", x@AIC)
    cat("\n\n")
  }
)

setMethod("show", "bdrm",
  function(object)
    getMethod("print", "bdrm")(object)
)

setMethod("print", "bdrm",
  function(x, dig = 3, ...) {
    ## title
    cat("Dose-response model\n\n")

    ## function call
    cat("Call:\n")
    print(x@call)

    ## model family & estimates
    family <-
      switch(x@family()$family,
             betapoisson = "Beta-Poisson",
             exponential = "Exponential",
             loglogistic = "Log-Logistic",
             logprobit = "Log-Probit",
             extremevalue = "Extreme value")
    cat("\n", family, " model coefficients:\n", sep = "")

    est <- summarize(x)
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

setGeneric("summarize",
  function(x, ...) {
    standardGeneric("summarize")
  }
)

setMethod("summarize", "ea",
  function(x, ...) {
    do.call(x@family()$summarize, as.list(x@coef))
  }
)

setMethod("summarize", "bea",
  function(x, ...) {
    posterior <- unlist(x@mcmc)
    return(data.frame(mean = mean(posterior),
                      sd   = sd(posterior),
                      "2.5%"  = quantile(posterior, .025),
                      "97.5%" = quantile(posterior, .975),
                      check.names = FALSE,
                      row.names = ""))
  }
)

setMethod("summarize", "bdrm",
  function(x, ...) {
    posterior <- as.matrix(x@mcmc)
    return(data.frame(mean = apply(posterior, 2, mean),
                      sd   = apply(posterior, 2, sd),
                      "2.5%"  = apply(posterior, 2, quantile, probs = 0.025),
                      "97.5%" = apply(posterior, 2, quantile, probs = 0.975),
                      check.names = FALSE,
                      row.names = x@family()$bayes$nodes))
  }
)

setGeneric("sim",
  function(x, n, ...) {
    standardGeneric("sim")
  }
)

setMethod("sim", "ea",
  function(x, n, ...) {
    do.call(x@family()$sim, c(as.list(x@coef), n))
  }
)

setMethod("sim", "bea",
  function(x, n, ...) {
    sample(unlist(x@mcmc), n, replace = TRUE)
  }
)

setMethod("sim", "drm",
  function(x, n, dose, ...) {
    do.call(x@family()$sim,
            list(n = n, mu = x@coef, Sigma = x@vcov, dose = dose))
  }
)

setMethod("sim", "bdrm",
  function(x, n, dose, ...) {
    pars <- apply(as.matrix(x@mcmc), 2, sample, size = n, replace = TRUE)
    do.call(x@family()$bayes$sim,
            list(pars = matrix(pars, ncol = 2),
                 dose = matrix(dose, ncol = 1)))
  }
)

setMethod("plot", "drm",
  function(x, y, se = TRUE, add = FALSE, n = NULL, min_log10dose = 0,
           xlab = "log10(dose)", ylab = "P(infection)",
           type = "l", lwd = 2, col = "red",
           se_pars = list(type = "l", lty = 2, lwd = 2, col = "blue"),
           sim_pars = list(type = "l", lty = 1, col = rgb(0, 0, 0, 0.1)),
           ...) {
    max_log10dose <- max(ceiling(log10(x@data$dose)))
    d <- seq(min_log10dose, max_log10dose, .1)
    p <- predict(x, 10^d)

    if (!add) {
      ## plot fitted dose-response curve
      plot(d, p[, 1],
           ylim = c(0, 1), xlim = c(min_log10dose, max_log10dose),
           xlab = xlab, ylab = ylab, 
           type = type, lwd = lwd, col = col, ...)

      ## plot observations
      points(log10(x@data$dose), x@data$x / x@data$n,
             cex = log(x@data$n))

    } else {
      ## plot fitted dose-response curve
      lines(d, p[, 1], type = type, lwd = lwd, col = col, ...)
    }

    ## plot SE curves
    if (se) {
      do.call(matlines,
              c(list(d, p[, 3:4]), se_pars))
    }

    ## plot simulated curves
    if (!is.null(n)) {
      do.call(matlines,
              c(list(d, sim(x, n = n, dose = 10^d)), sim_pars))
    }
  }
)

setMethod("plot", "bdrm",
  function(x, y, n = 100, add = FALSE, ...) {
    max_log10dose <- max(ceiling(log10(x@data$dose)))
    d <- seq(0, max_log10dose, .1)
    p <- sim(x, n = n, 10^d)

    if (!add) {
      ## plot fitted dose-response curves
      matplot(d, p,
              xlab = "log10(dose)", ylab = "P(infection)", 
              ylim = c(0, 1), xlim = c(0, max_log10dose),
              col = rgb(0, 0, 0, .05), lty = 1, type = "l", ...)

      ## plot mean dose-response curve
      lines(d, apply(p, 1, mean), lwd = 2)

      ## plot observations
      points(log10(x@data$dose), x@data$x / x@data$n,
            cex = log(x@data$n))
    } else {
      ## plot fitted dose-response curve
      matlines(d, p, lty = 1, type = "l", ...)
    }
  }
)

setMethod("predict", "drm",
  function(object, dose, conf_level = 0.95, ...) {
    pred <- do.call(object@family()$predict,
                    list(object@coef, object@vcov, dose))
    pred <- matrix(pred, ncol = 2)
    pred_ci <- logit_transformed_ci(est = pred[, 1], se = pred[, 2],
                                    conf_level = conf_level)
    out <- cbind(pred, pred_ci)
    colnames(out) <- c("estimate", "se",
                       paste0(100 * (1 - conf_level) / 2, "%"),
                       paste0(100 * (1 - (1 - conf_level) / 2), "%"))
    rownames(out) <- paste0("dose", dose)
    return(out)
  }
)
