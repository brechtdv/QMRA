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