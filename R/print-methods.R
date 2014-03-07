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
  cat("Exposure assessment from", from, "data\n\n")

  ## model estimates
  cat("Individual model estimates:\n")
  print(x$mod)

  ## model estimates
  cat("\nModel average:\n")
  print(x$avg)
}