## Print JAGS model
print.JAGS_model <-
function(x, ...){
  l <- length(x)
  spacer <- 0
  for (i in seq(l)){
    if (substr(x[i], nchar(x[i]), nchar(x[i])) == "}")
      spacer <- spacer - 1
    cat(rep(" ", 2 * spacer), x[i], "\n", sep = "")
    if (substr(x[i], nchar(x[i]), nchar(x[i])) == "{")
      spacer <- spacer + 1
  }
}

## Print model average
print.avg <-
function(x, ...) {
  ## data type
  from <-
    c("count", "concentration", "presence/absence")[
      which(x$fun == c("ea_count", "ea_conc", "ea_presence"))]
  if (length(from) != 0) {
    cat("Exposure assessment from", from, "data\n\n")
  } else {
    cat("Dose-response assessment\n\n")
  }

  ## model estimates
  cat("Individual model estimates:\n")
  print(x$mod)

  ## model estimates
  cat("\nModel average:\n")
  print(x$avg)
}

## Print Goodness of Fit
print.gof <-
function(x, ...) {
  if (is.null(x)) {
    cat("Goodness of fit statistics not yet implemented for this function.\n")

  } else {
    cat("Chi-square goodness of fit test\n\n")
    cat("Saturated likelihood:",
        x$L0[1], "on", x$L0[2], "degree(s) of freedom\n")
    cat("    Model likelihood:",
        x$L1[1], "on", x$L1[2], "degree(s) of freedom\n\n")

    df <- x$L0[2] - x$L1[2]
    chisq <- -2 * log(x$L1[1] / x$L0[1])
    p     <- pchisq(chisq, df, lower.tail = FALSE)

    out <- data.frame(chisq, df, p)
    rownames(out) <- ""
    colnames(out) <- c("ChiSq", "df", "Pr>ChiSq")
    print(out)
  }
}