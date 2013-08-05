## EXPOSURE ASSESSMENT

## COUNT DATA --------------------------------------------------------------
ea_count <-
function(x, q = 1, data,
         family = c("poisson", "negbin", "poislognorm", "poisinvgauss", "poisgeninvgauss")){

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

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family()$minloglik,
        start = family()$start,
        fixed = list(x = x, q = q))

  ## AIC = -2*logLik + 2*npar
  AIC <- summary(MLE)@m2logL + 2 * family()$npar

  ## Likelihood Ratio test
  LL0 <- -logLik(MLE)
  LL1 <- poisson()$minloglik(x / q, x = x, q = q)
  LRT <- 2 * (LL0 - LL1)
  p <- pchisq(LRT, length(x) - family()$npar, lower.tail = FALSE)
  lrt <- list("LL0" = LL0, "LL1" = LL1, "LRT" = LRT, "p-chisq" = p)

  ## create 'ea' object
  ea_fit <-
    new("ea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        estimates = lrt,
        AIC = AIC,
        mle = MLE)

  ## return 'ea' object
  return(ea_fit)
}

## CONCENTRATION DATA ------------------------------------------------------
ea_conc <-
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

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family(x, d)$minloglik,
        start = family(x, d)$start,
        fixed = list(x = x, d = d))

  ## AIC = -2*logLik + 2*npar
  AIC <- summary(MLE)@m2logL + 2 * family(x, d)$npar

  ## create 'ea' object
  ea_fit <-
    new("ea",
        call = call,
        data = data.frame(x = x, d = d),
        family = family,
        estimates = list(),
        AIC = AIC,
        mle = MLE)

  ## return 'ea' object
  return(ea_fit)
}

## PRESENCE/ABSENCE DATA ---------------------------------------------------
ea_presence <-
function(x, q = 1, replicates = rep(1, length(x)), data,
         family = c("poisson")){

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

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family()$minloglik_bernoulli,
        start = family()$start,
        fixed = list(x = x, q = q))

  ## AIC = -2*logLik + 2*npar
  AIC <- summary(MLE)@m2logL + 2 * family()$npar

  ## create 'ea' object
  ea_fit <-
    new("ea",
        call = call,
        data = data.frame(x = x, q = q),
        family = family,
        estimates = list(),
        AIC = AIC,
        mle = MLE)

  ## return 'ea' object
  return(ea_fit)
}