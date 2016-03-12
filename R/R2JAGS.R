R2JAGS <-
function(model, data, nchains, inits, burnin, nodes, update, verbose) {
  ## disable JAGS progress bars
  if (!verbose) {  
    old.pb <- options("jags.pb")
    on.exit(options(old.pb)) 
    options(jags.pb = "none")
  }

  ## format model
  wrap(model)

  ## define & initialize model
  mod <-
    jags.model(file = "modelTempFile.txt",
               data = data,
               inits = inits,
               n.chains = nchains,
               n.adapt = 1000,
               quiet = !verbose)

  ## delete 'modelTempFile.txt'
  unlink("modelTempFile.txt")

  ## burn-in
  update(mod, n.iter = burnin)

  ## samples
  samples <- coda.samples(mod, nodes, n.iter = update, thin = 1)

  ## Deviance
  dic <-
    dic.samples(mod, n.iter = update, thin = 1, type = "pD")

  ## Return results
  return(list(mcmc.list = samples, dic = dic))
}
