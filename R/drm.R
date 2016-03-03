## DOSE-RESPONSE MODELLING

drm <-
function(x, n, dose, data,
         model = c("betapoisson", "bp", 
                   "exponential", "exp", 
                   "loglogistic", "ll",
                   "logprobit", "lp",
                   "extremevalue", "ev"),
         start = NULL, ...) {

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

  ## define optimizer start values
  ## use default value if no value defined by user
  if (is.null(start)) {
    start <- family()$start
  }

  ## get maximum likelihood estimate
  MLE <-
    mle(minuslogl = family()$minloglik,
        start = start,
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
