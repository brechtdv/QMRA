## MODEL AVERAGING

## MAIN FUNCTION --------------------------------------------------------------
avg <-
function(...) {
  ## group models
  models <- list(...)

  ## check if all elements are same class...

  ## derive generating function
  fun <- models[[1]]@call[[1]]

  ## derive weights
  aic <- sapply(models, function(x) x@AIC)
  delta <- aic - min(aic)
  wghts <- exp(-.5 * delta) / sum(exp(-.5 * delta))

  ## derive mean and variance
  m <- sapply(models, function(x) summarize(x)$mean)
  v <- sapply(models, function(x) summarize(x)$sd) ^ 2

  ## calculate average mean and variance
  m_avg <- sum(wghts * m)
  s_avg <- sum(wghts * sqrt(v + (m - m_avg)^2))

  ## calculate lower and upper bound ~ lognormal approximation
  C <- exp(qnorm(.975) * sqrt(log(1 + (s_avg / m_avg)^2)))
  cnfi <- c(m_avg / C, m_avg * C)

  ## define function output
  out <-
    list(fun = fun,
         mod = data.frame(mean = m, sd = sqrt(v), AIC = aic, weight = wghts),
         avg = data.frame(mean = m_avg, sd = s_avg,
                          "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                          check.names = FALSE, row.names = ""))
  class(out) <- "avg"

  ## return output
  return(out)
}

## Bayesian model averaging ------------------------------------------------
avg_bea <-
function(...) {
  ## group models
  models <- list(...)

  ## check if all elements are same class...

  ## derive generating function
  fun <- models[[1]]@call[[1]]

  ## derive weights
  dic <- sapply(models, function(x) sum(unlist(x@diagnostics$DIC[-3])))
  delta <- dic - min(dic)
  wghts <- exp(-.5 * delta) / sum(exp(-.5 * delta))

  ## derive mean and variance
  m <- sapply(models, function(x) summarize(x)$mean)
  v <- sapply(models, function(x) summarize(x)$sd) ^ 2

  ## calculate average mean and variance
  m_avg <- sum(wghts * m)
  s_avg <- sum(wghts * sqrt(v + (m - m_avg)^2))

  ## calculate lower and upper bound ~ lognormal approximation
  C <- exp(qnorm(.975) * sqrt(log(1 + (s_avg / m_avg)^2)))
  cnfi <- c(m_avg / C, m_avg * C)

  ## define function output
  out <-
    list(fun = fun,
         mod = data.frame(mean = m, sd = sqrt(v), DIC = dic, weight = wghts),
         avg = data.frame(mean = m_avg, sd = s_avg,
                          "2.5%" = cnfi[1], "97.5%" = cnfi[2],
                          check.names = FALSE, row.names = ""))
  class(out) <- "avg"

  ## return output
  return(out)
}