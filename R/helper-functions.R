logit <-
function(x) {
  log(x / (1 - x))
}

expit <-
function(x) {
  exp(x) / (1 + exp(x))
}

delta_method <-
function(g, coef, vcov) {
  n <- length(coef)
  syms <- paste("x", seq(n), sep = "")
  for (i in seq(n)) assign(syms[i], coef[i])
  gdashmu <-
    t(sapply(g, function(form) {
        as.numeric(attr(eval(deriv(form, syms)), "gradient"))
    }))
  new_covar <- gdashmu %*% vcov %*% t(gdashmu)
  new_se <- sqrt(diag(new_covar))
  return(new_se)
}

logit_transformed_ci <-
function(est, se, conf_level) {
  eta <- log(est / (1 - est))
  se <- se / (est * (1 - est))
  ll <- eta + qnorm((1 - conf_level) / 2) * se
  ul <- eta - qnorm((1 - conf_level) / 2) * se
  lwr <- exp(ll) / (1 + exp(ll))
  upr <- exp(ul) / (1 + exp(ul))
  return(cbind(lwr, upr))
}
