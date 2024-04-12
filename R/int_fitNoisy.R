###   Functions to fit the Gaussian noise  model but slightly improved model fitting algorithm
# rating criteria are enforced to be ordered, but there can be an overlap between rating criteria and type 1 criterion
# use the best fitting set of parameters of the initial grid search
# include arguenmts to set how many starting values should be tried and how often the fitting routine should be restarted

fitNoisy <-
  function(N_SA_RA, N_SA_RB, N_SB_RA, N_SB_RB,
           nInits, nRestart, nRatings, nCond, nTrials){

    # coarse grid search to find promising initial values

    temp <- expand.grid(maxD =  seq(1, 5, 1),
                        theta = seq(-1/2,1/2, 1/2),
                        tauMin =  c(.1, .3, 1),  # position of the most conservative confidence criteria with respect to theta
                        tauRange = seq(1, 5, 1),  #  position of the most liberal confidence criterion with respect to theta
                        sigma = c(.1, .3, 1, 3)) # noise parameter

    # number of parameters:
    # nCond -1 sensitivity parameters
    # 1 type 1 bias parameter
    # 2 * (nRatings -1) confidence criteria
    # 1 confidence noise parameters

    inits <- data.frame(matrix(data=NA, nrow= nrow(temp), ncol = nCond + nRatings*2))
    if(nCond==1)  {
      inits[,1] <-  log(temp$maxD)
    } else{
      inits[,1:(nCond)] <-  log(t(mapply(function(maxD) diff(seq(0, maxD, length.out = nCond+1)), temp$maxD)))
    }
    if (nRatings > 3){
      inits[,(nCond+1):(nCond+nRatings-2)] <-
        log(t(mapply(function(tauMin, tauRange) diff(seq(-tauRange-tauMin, -tauMin, length.out=nRatings-1)),
                     temp$tauMin, temp$tauRange)))
      inits[,(nCond+nRatings+2):(nCond + nRatings*2-1)] <-
        log(t(mapply(function(tauMin, tauRange) diff(seq(tauMin, tauMin+tauRange, length.out=nRatings-1)),
                     temp$tauMin, temp$tauRange)))
    }
    if (nRatings == 3){
      inits[,(nCond+1):(nCond+nRatings-2)] <-
        log(mapply(function(tauMin, tauRange) diff(seq(-tauRange-tauMin, -tauMin, length.out=nRatings-1)),
                     temp$tauMin, temp$tauRange))
      inits[,(nCond+nRatings+2):(nCond + nRatings*2-1)] <-
        log(mapply(function(tauMin, tauRange) diff(seq(tauMin, tauMin+tauRange, length.out=nRatings-1)),
                     temp$tauMin, temp$tauRange))
    }
    inits[,nCond+(nRatings-1)] <- -temp$tauMin
    inits[,nCond+nRatings] <- temp$theta
    inits[,nCond+(nRatings+1)] <- temp$tauMin
    inits[,(nCond + nRatings*2)] <- log(temp$sigma)

    logL <- apply(inits, MARGIN = 1,
                  function(p) try(ll_Noisy(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond), silent = TRUE))
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    inits <- inits[1:nInits,]

    noFitYet <- TRUE
    for (i in 1:nInits){
      m <- try(optim(par =  inits[i,],
                     fn = ll_Noisy, gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                     control = list(maxit = 10^4, reltol = 10^-4)))

      if (is.list(m)){
        for(j in 2:nRestart){
          try(m <- optim(par = m$par,
                         fn = ll_Noisy, gr = NULL,
                         N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                         N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                         control = list(maxit = 10^6, reltol = 10^-8)))

        }
        if (noFitYet) {
          fit <- m
          noFitYet <- FALSE
        } else {
          if (m$value < fit$value) fit <- m
        }
      }
    }

    res <-  data.frame(matrix(nrow=1, ncol=0))
    if(!inherits(fit, "try-error")){

      k <- length(fit$par)
      res[paste("d_",1:nCond, sep="")] <-  as.vector(cumsum(exp(fit$par[1:(nCond)])))
      res$c <-  as.vector(fit$par[nCond+nRatings])
      res[,paste("theta_minus.",(nRatings-1):1, sep="")] <-
        c(as.vector(fit$par[nCond+nRatings-1] - rev(cumsum(c(exp(fit$par[(nCond+1):(nCond+nRatings-2)]))))),
          as.vector(fit$par[nCond+nRatings-1]))
      res[,paste("theta_plus.",1:(nRatings-1), sep="")] <-
        c(as.vector(fit$par[nCond+nRatings+1]),
          as.vector(fit$par[nCond+nRatings+1]) +
            as.vector(cumsum(c(exp(fit$par[(nCond+nRatings+2):(nCond + nRatings*2-1)])))))


      res$sigma <- exp(fit$par[nCond + nRatings*2])

      res$negLogLik <- fit$value
      res$N <- nTrials
      res$k <- k
      res$BIC <-  2 * fit$value + k * log(nTrials)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(nTrials-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }
