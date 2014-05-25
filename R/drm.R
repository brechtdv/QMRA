## DOSE-RESPONSE MODELLING

drm <-
function(x, n, dose, data,
         model = c("bp", "betapoisson", "exp", "exponential",
                   "ll", "loglogistic", "lp", "logprobit",
                   "ev", "extremevalue"), ...){

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

  ## check model
  model <- match.arg(model)
  if (any(model == c("bp", "exp", "ll", "lp", "ev"))) {
    model <-
      c("betapoisson", "exponential",
        "loglogistic", "logprobit", "extremevalue")[
        which(model == c("bp", "exp", "ll", "lp", "ev"))]
  }

  ## obtain 'family()' function
  family <- get(model)
  #family <- getFromNamespace(model, "QMRA")

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family()$minloglik,
        start = family()$start,
        fixed = list(x = x, n = n, d = d), ...)

  ## AIC = -2*logLik + 2*npar
  AIC <- summary(MLE)@m2logL + 2 * family()$npar

  ## Chi-square goodness of fit
  p0 <- x / n
  L0 <- prod(p0 ^ x * (1 - p0) ^ (n - x))
  L1 <- exp(summary(MLE)@m2logL / -2)
  gof <- list(L0 = c(L0, length(x)), L1 = c(L1, family()$npar))
  class(gof) <- "gof"

  ## create 'drm' object
  drm_fit <-
    new("drm",
        call = call,
        data = data.frame(x = x, n = n, dose = d),
        family = family,
        AIC = AIC,
        gof = gof,
        mle = MLE)

  ## return 'drm' object
  return(drm_fit)
}
