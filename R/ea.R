## EXPOSURE ASSESSMENT

## COUNT DATA --------------------------------------------------------------
ea_count <-
function(x, q = 1, data,
         model = c("poisson", "p",
                   "negbin", "nb",
                   "poislognorm", "pln",
                   "poisinvgauss", "pig",
                   "poisgeninvgauss", "pgig"), ...){

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

  ## check 'model'
  model <- match.arg(model)
  full <- c("poisson", "negbin", "poislognorm",
            "poisinvgauss", "poisgeninvgauss")
  short <- c("p", "nb", "pln", "pig", "pgig")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## fit poisson-lognormal as generalized linear mixed model
  if (model == "poislognorm") {
    id <- seq_along(x)
    fit_glmm <-
      glmer(x ~ 1 + (1 | id), offset = log(q), family = stats::poisson, ...)

    MLE <- new("mle",
               coef = c(mu_log = fit_glmm@beta, sd_log = fit_glmm@theta))
    AIC <- AIC(fit_glmm)
    gof <- list()
    class(gof) <- "gof"

  ## fit other models via likelihood
  } else {
    ## get maximum likelihood estimate
    MLE <-
      mle(minuslogl = family()$minloglik,
          start = family()$start,
          fixed = list(x = x, q = q), ...)

    ## AIC = -2*logLik + 2*npar
    AIC <- summary(MLE)@m2logL + 2 * family()$npar

    ## Chi-square goodness of fit
    L0 <- exp(-poisson()$minloglik(x / q, x = x, q = q))
    L1 <- exp(logLik(MLE))
    gof <- list(L0 = c(L0, length(x)), L1 = c(L1, family()$npar))
    class(gof) <- "gof"
  }

  ## create 'ea' object
  ea_fit <-
    new("ea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        gof = gof,
        AIC = AIC,
        mle = MLE)

  ## return 'ea' object
  return(ea_fit)
}

## CONCENTRATION DATA ------------------------------------------------------
ea_conc <-
function(x, d, data,
         model = c("gamma", "g",
                   "lognorm", "ln",
                   "weibull", "w",
                   "invgauss", "ig"), ...){

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

  ## check 'model'
  model <- match.arg(model)
  full <- c("gamma", "lognorm", "weibull", "invgauss")
  short <- c("g", "ln", "w", "ig")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family(x, d)$minloglik,
        start = family(x, d)$start,
        fixed = list(x = x, d = d), ...)

  ## AIC = -2*logLik + 2*npar
  AIC <- summary(MLE)@m2logL + 2 * family(x, d)$npar

  ## Chi-square goodness of fit
  ## How to implement???
  gof <- NULL
  class(gof) <- "gof"

  ## create 'ea' object
  ea_fit <-
    new("ea",
        call = call,
        data = data.frame(x = x, d = d),
        family = family,
        gof = gof,
        AIC = AIC,
        mle = MLE)

  ## return 'ea' object
  return(ea_fit)
}

## PRESENCE/ABSENCE DATA ---------------------------------------------------
ea_presence <-
function(x, q = 1, replicates = rep(1, length(x)), data,
         model = c("poisson", "p"), ...){

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

  ## check 'model'
  model <- match.arg(model)
  full <- c("poisson")
  short <- c("p")
  if (model %in% short) model <- full[match(model, short)]

  ## obtain 'family()' function
  family <- getFromNamespace(model, "QMRA")

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family()$minloglik_bernoulli,
        start = family()$start,
        fixed = list(x = x, q = q), ...)

  ## AIC = -2*logLik + 2*npar
  AIC <- summary(MLE)@m2logL + 2 * family()$npar

  ## Chi-square goodness of fit
  L0 <- exp(-poisson()$minloglik_bernoulli(x, x, q))
  L1 <- exp(logLik(MLE))
  gof <- list(L0 = c(L0, length(x)), L1 = c(L1, family()$npar))
  class(gof) <- "gof"

  ## create 'ea' object
  ea_fit <-
    new("ea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        gof = gof,
        AIC = AIC,
        mle = MLE)

  ## return 'ea' object
  return(ea_fit)
}