###   Functions to fit the standard SDT model
# the model does not assume any parameter to measure metacognition
# rating criteria are enforced to be ordered
# sensitivity is / are forced to be ordered and to be greater than 0

fitSDT <-
  function(N_SA_RA, N_SA_RB, N_SB_RA, N_SB_RB,
           nInits, nRestart, nRatings, nCond){

    # coarse grid search to find promising initial values

    temp <- expand.grid(maxD =  seq(1,5,1),
                        theta = seq(-1/2,1/2, 1/2),
                        tauMin =  seq(.1,1, length.out=4),  # distance of the first rating criterion from the task crtierion for stimulus A
                        tauRange = seq(1,5,length.out=5)) # range of rating criteria stimulus B

    inits <- data.frame(matrix(data=NA, nrow= nrow(temp), ncol = nCond + nRatings*2 - 1))
    if(nCond==1)  {
      inits[,1] <-  log(temp$maxD)  }
    else{
      inits[,1:(nCond)] <-  log(t(mapply(function(maxD) diff(seq(0, maxD, length.out = nCond+1)), temp$maxD)))
    }
    if (nRatings > 3){
      inits[,(nCond+1):(nCond+nRatings-2)] <-
        log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                     temp$tauRange)))
      inits[,(nCond+nRatings+2):ncol(inits)] <-
        log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                     temp$tauRange)))
    }
    if (nRatings == 3){
      inits[,(nCond+1):(nCond+nRatings-2)] <-
        log(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                   temp$tauRange))
      inits[,(nCond+nRatings+2):ncol(inits)] <-
        log(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                   temp$tauRange))
    }
    inits[,nCond+(nRatings-1)] <- log(temp$tauMin)
    inits[,nCond+nRatings] <- temp$theta
    inits[,nCond+(nRatings+1)] <- log(temp$tauMin)

    logL <- apply(inits, MARGIN = 1,
                  function(p) try(llSDT(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond), silent = TRUE))
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    inits <- inits[1:nInits,]

    noFitYet <- TRUE
    for (i in 1:nInits){
      m <- try(optim(par =  inits[i,],
                     fn = llSDT, gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                     control = list(maxit = 10^4, reltol = 10^-4)))

      if (is.list(m)){
        for(j in 2:nRestart){
          try(m <- optim(par = m$par,
                         fn = llSDT, gr = NULL,
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
      N <- length(ratings)

      res[paste("d_",1:nCond, sep="")] <-  as.vector(cumsum(exp(fit$par[1:(nCond)])))
      res$c <-  as.vector(fit$par[nCond+nRatings])
      res[,paste("theta_minus.",(nRatings-1):1, sep="")] <-
        as.vector(fit$par[nCond+nRatings] - rev(cumsum(c(exp( fit$par[(nCond+1):(nCond+nRatings-1)])))))
      res[,paste("theta_plus.",1:(nRatings-1), sep="")] <-
        as.vector(fit$par[nCond+nRatings] + cumsum(c(exp(fit$par[(nCond+nRatings+1):(nCond + nRatings*2-1)]))))


      res$negLogLik <- fit$value
      res$N <- N
      res$k <- k
      res$BIC <-  2 * fit$value + k * log(N)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }

