# modified version of the BESTmcmc function in the package BEST
### modified to return MCMC chains for posterior processing in analyses that pool
##### results from imputed datasets

BESTmcmc_modified <- 
function (y1, y2 = NULL, priors = NULL, doPriorsOnly = FALSE, 
          numSavedSteps = 1e+05, thinSteps = 1, burnInSteps = 1000, 
          verbose = TRUE, rnd.seed = NULL, parallel = NULL) 
{
  if (doPriorsOnly && verbose) 
    cat("Warning: The output shows the prior distributions,\n      NOT the posterior distributions for your data.\n")
  nCores <- detectCores()
  if (!is.null(parallel) && parallel && nCores < 4) {
    if (verbose) 
      warning("Not enough cores for parallel processing, running chains sequentially.")
    parallel <- FALSE
  }
  if (is.null(parallel)) 
    parallel <- nCores > 3
  y <- c(y1, y2)
  if (!all(is.finite(y))) 
    stop("The input data include NA or Inf.")
  if (length(unique(y)) < 2 && (is.null(priors) || is.null(priors$muSD) || 
                                is.null(priors$sigmaMode) || is.null(priors$sigmaSD))) 
    stop("If priors are not specified, data must include at least 2 (non-equal) values.")
  if (!is.null(priors)) {
    if (!is.list(priors)) {
      if (is.numeric(priors)) {
        stop("'priors' is now the 3rd argument; it must be a list (or NULL).")
      }
      else {
        stop("'priors' must be a list (or NULL).")
      }
    }
    nameOK <- names(priors) %in% c("muM", "muSD", "sigmaMode", 
                                   "sigmaSD", "nuMean", "nuSD")
    if (!all(nameOK)) 
      stop("Invalid items in prior specification: ", paste(sQuote(names(priors)[!nameOK]), 
                                                           collapse = ", "))
    if (!all(sapply(priors, is.numeric))) 
      stop("All items in 'priors' must be numeric.")
    if (!is.null(priors$muSD) && priors$muSD <= 0) 
      stop("muSD must be > 0")
  }
  if (is.null(rnd.seed)) 
    rnd.seed <- floor(runif(1, 1, 10000))
  if (is.null(priors)) {
    dataForJAGS <- list(muM = mean(y), muP = 1e-06 * 1/sd(y)^2, 
                        sigmaLow = sd(y)/1000, sigmaHigh = sd(y) * 1000)
  }
  else {
    priors0 <- list(muM = mean(y), muSD = sd(y) * 5, sigmaMode = sd(y), 
                    sigmaSD = sd(y) * 5, nuMean = 30, nuSD = 30)
    priors0 <- modifyList(priors0, priors)
    sigmaShRa <- gammaShRaFromModeSD(mode = priors0$sigmaMode, 
                                     sd = priors0$sigmaSD)
    nuShRa <- gammaShRaFromMeanSD(mean = priors0$nuMean, 
                                  sd = priors0$nuSD)
    dataForJAGS <- list(muM = priors0$muM, muP = 1/priors0$muSD^2, 
                        Sh = sigmaShRa$shape, Ra = sigmaShRa$rate)
    if (!is.null(y2)) {
      fixPrior <- function(x) {
        if (length(x) < 2) 
          x <- rep(x, 2)
        return(x)
      }
      dataForJAGS <- lapply(dataForJAGS, fixPrior)
    }
    dataForJAGS$ShNu <- nuShRa$shape
    dataForJAGS$RaNu <- nuShRa$rate
  }
  modelFile <- file.path(tempdir(), "BESTmodel.txt")
  if (is.null(priors)) {
    if (is.null(y2)) {
      modelString = "\n      model {\n        for ( i in 1:Ntotal ) {\n          y[i] ~ dt( mu , tau , nu )\n        }\n        mu ~ dnorm( muM , muP )\n        tau <- 1/pow( sigma , 2 )\n        sigma ~ dunif( sigmaLow , sigmaHigh )\n        nu <- nuMinusOne+1\n        nuMinusOne ~ dexp(1/29)\n      }\n      "
    }
    else {
      modelString <- "\n      model {\n        for ( i in 1:Ntotal ) {\n          y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )\n        }\n        for ( j in 1:2 ) {\n          mu[j] ~ dnorm( muM , muP )\n          tau[j] <- 1/pow( sigma[j] , 2 )\n          sigma[j] ~ dunif( sigmaLow , sigmaHigh )\n        }\n        nu <- nuMinusOne+1\n        nuMinusOne ~ dexp(1/29)\n      }\n      "
    }
  }
  else {
    if (is.null(y2)) {
      modelString = "\n      model {\n        for ( i in 1:Ntotal ) {\n          y[i] ~ dt( mu , tau , nu )\n        }\n        mu ~ dnorm( muM[1] , muP[1] )\n        tau <- 1/pow( sigma , 2 )\n         sigma ~ dgamma( Sh[1] , Ra[1] )\n        nu ~ dgamma( ShNu , RaNu ) # prior for nu\n      }\n      "
    }
    else {
      modelString <- "\n      model {\n        for ( i in 1:Ntotal ) {\n          y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )\n        }\n        for ( j in 1:2 ) {\n          mu[j] ~ dnorm( muM[j] , muP[j] )\n          tau[j] <- 1/pow( sigma[j] , 2 )\n          sigma[j] ~ dgamma( Sh[j] , Ra[j] )\n        }\n        nu ~ dgamma( ShNu , RaNu ) # prior for nu\n      }\n      "
    }
  }
  writeLines(modelString, con = modelFile)
  if (!doPriorsOnly) 
    dataForJAGS$y <- y
  dataForJAGS$Ntotal <- length(y)
  if (!is.null(y2)) 
    dataForJAGS$x <- c(rep(1, length(y1)), rep(2, length(y2)))
  if (is.null(y2)) {
    mu = mean(y1)
    sigma = sd(y1)
  }
  else {
    mu = c(mean(y1), mean(y2))
    sigma = c(sd(y1), sd(y2))
  }
  initsList0 <- list(mu = mu, sigma = sigma, .RNG.seed = rnd.seed)
  if (is.null(priors)) {
    initsList0$nuMinusOne <- 4
  }
  else {
    initsList0$nu <- 5
  }
  initsList <- list(c(initsList0, .RNG.name = "base::Wichmann-Hill"), 
                    c(initsList0, .RNG.name = "base::Marsaglia-Multicarry"), 
                    c(initsList0, .RNG.name = "base::Super-Duper"))
  codaSamples <- jags.basic(data = dataForJAGS, inits = initsList, 
                            parameters.to.save = c("mu", "sigma", "nu"), model.file = modelFile, 
                            n.chains = 3, n.adapt = 500, n.iter = ceiling((numSavedSteps * 
                                                                             thinSteps)/3 + burnInSteps), n.burnin = burnInSteps, 
                            n.thin = thinSteps, modules = NULL, parallel = parallel, 
                            DIC = FALSE, seed = rnd.seed, save.model=TRUE,verbose = verbose)
  
  return(codaSamples)
}
environment(BESTmcmc_modified) <- asNamespace('BEST')
