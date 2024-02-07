###  Fit the independent truncated Gaussian model
### there exist two versions:
### - ITGcm is identical to the version of meta-d' proposed by Maniscalco and Lau (2012)
### - ITGc is identical to HMetad proposed by Fleming (2017)


fitITGcm <-
  function(ratings, stimulus, correct, condition, nInits = 5, nRestart = 4){

    A <- levels(stimulus)[1]
    B <- levels(stimulus)[2]
    nRatings <- length(levels(ratings))
    nCond <- length(levels(condition))

    N_SA_RA <- table(condition[stimulus == A & correct == 1],
                     ratings[stimulus == A & correct == 1])[,nRatings:1] + 0.001
    N_SA_RB <- table(condition[stimulus == A & correct == 0],
                     ratings[stimulus == A & correct == 0]) + 0.001
    N_SB_RA <- table(condition[stimulus == B & correct == 0],
                     ratings[stimulus == B & correct == 0])[,nRatings:1] + 0.001
    N_SB_RB <- table(condition[stimulus == B & correct == 1],
                     ratings[stimulus == B & correct == 1]) + 0.001

    # search for initial values using a coarse search grind
    temp <- expand.grid(maxD =  seq(1, 5, 1),
                        theta = seq(-1/2,1/2, 1/2),
                        tauMin =  c(.1, .3, 1),  # position of the most conservative confidence criteria with respect to theta
                        tauRange = seq(1, 5, 1),  #  position of the most liberal confidence criterion with respect to theta
                        m = c(.1, .3, 1, 3)) # metacognitive efficiency

    # number of parameters:
    # nCond -1 sensitivity parameters
    # 1 type 1 bias parameter
    # 2 * (nRatings - 1) confidence criteria
    # m-ratio parameter
    #
    inits <- data.frame(matrix(data=NA, nrow= nrow(temp),
                               ncol = nCond + nRatings*2 ))
    if(nCond==1)  {
      inits[,1] <-  log(temp$maxD)  }
    else{
      inits[,1:(nCond)] <-  log(t(mapply(function(maxD) diff(seq(0, maxD, length.out = nCond+1)), temp$maxD)))
    }
    inits[,(nCond+1):(nCond+nRatings-2)] <-
      log(t(mapply(function(tauMin, tauRange) diff(seq(-tauRange-tauMin, -tauMin, length.out=nRatings-1)),
                   temp$tauMin, temp$tauRange)))
    inits[,nCond+(nRatings-1)] <- log(temp$tauMin)
    inits[,nCond+nRatings] <- temp$theta
    inits[,nCond+(nRatings+1)] <- log(temp$tauMin)
    inits[,(nCond+nRatings+2):(nCond + nRatings*2-1)] <-
      log(t(mapply(function(tauMin, tauRange) diff(seq(tauMin, tauMin+tauRange, length.out=nRatings-1)),
                   temp$tauMin, temp$tauRange)))

    inits[,(nCond + nRatings*2)] <- log(temp$m)

    logL <- apply(inits, MARGIN = 1,
                  function(p) try(ll_Mratio(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond), silent = TRUE))
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    inits <- inits[1:nInits,]

    noFitYet <- TRUE
    for (i in 1:nInits){

      m <- try(optim(par =  inits[i,],
                     fn = ll_Mratio, gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                     control = list(maxit = 10^4, reltol = 10^-4)))

      if (is.list(m)){
        for(j in 2:nRestart){
          try(m <- optim(par = m$par,
                         fn = ll_Mratio, gr = NULL,
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
        # c(as.vector(fit$par[nCond+nRatings-1] - rev(cumsum(c(exp(fit$par[(nCond+1):(nCond+nRatings-2)]))))),  as.vector(fit$par[nCond+nRatings-1]))
        exp(fit$par[nCond + nRatings*2]) * as.vector(fit$par[nCond+nRatings]) -
        rev( cumsum(c(exp(fit$par[(nCond+1):(nCond+nRatings-1)]))))

      res[,paste("theta_plus.",1:(nRatings-1), sep="")] <-
        # c(as.vector(fit$par[nCond+nRatings+1]), as.vector(fit$par[nCond+nRatings+1]) + as.vector(cumsum(c(exp(fit$par[(nCond+nRatings+2):(nCond + nRatings*2-1)])))))
        exp(fit$par[nCond + nRatings*2]) * as.vector(fit$par[nCond+nRatings]) +
        cumsum(c(exp(fit$par[(nCond+nRatings+1):(nCond + nRatings*2-1)])))

      res$m <- exp(fit$par[nCond + nRatings*2])

      res$negLogLik <- fit$value
      res$N <- N
      res$k <- k
      res$BIC <-  2 * fit$value + k * log(N)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }

fitITGc <-
  function(ratings, stimulus, correct, condition, nInits = 5, nRestart = 4){
    if(!is.factor(condition)) stop ("condition should be a factor!")
    if(!is.factor(ratings)) stop ("ratings should be a factor!")
    if(!is.factor(stimulus )|| length(levels(stimulus)) != 2) {
      stop("stimulus should be a factor with 2 levels")
    }
    if(!all(correct %in% c(0,1))) stop("correct should be 1 or 0")

    A <- levels(stimulus)[1]
    B <- levels(stimulus)[2]
    nRatings <- length(levels(ratings))
    nCond <- length(levels(condition))

    N_SA_RA <- table(condition[stimulus == A & correct == 1],
                     ratings[stimulus == A & correct == 1])[,nRatings:1] + 0.001
    N_SA_RB <- table(condition[stimulus == A & correct == 0],
                     ratings[stimulus == A & correct == 0]) + 0.001
    N_SB_RA <- table(condition[stimulus == B & correct == 0],
                     ratings[stimulus == B & correct == 0])[,nRatings:1] + 0.001
    N_SB_RB <- table(condition[stimulus == B & correct == 1],
                     ratings[stimulus == B & correct == 1]) + 0.001

    # search for inital values using a coarse search grind
    temp <- expand.grid(maxD =  seq(1, 5, 1),
                        theta = seq(-1/2,1/2, 1/2),
                        tauMin =  c(.1, .3, 1),  # position of the most conservative confidence criteria with respect to theta
                        tauRange = seq(1, 5, 1),  #  position of the most liberal confidence criterion with respect to theta
                        m = c(.1, .3, 1, 3)) # metacognitive efficiency

    inits <- data.frame(matrix(data=NA, nrow= nrow(temp),
                               ncol = nCond + nRatings*2 ))
    if(nCond==1)  {
      inits[,1] <-  log(temp$maxD)  }
    else{
      inits[,1:(nCond)] <-  log(t(mapply(function(maxD) diff(seq(0, maxD, length.out = nCond+1)), temp$maxD)))
    }
    inits[,(nCond+1):(nCond+nRatings-2)] <-
      log(t(mapply(function(tauMin, tauRange) diff(seq(-tauRange-tauMin, -tauMin, length.out=nRatings-1)),
                   temp$tauMin, temp$tauRange)))
    inits[,nCond+(nRatings-1)] <- log(temp$tauMin)
    inits[,nCond+nRatings] <- temp$theta
    inits[,nCond+(nRatings+1)] <- log(temp$tauMin)
    inits[,(nCond+nRatings+2):(nCond + nRatings*2-1)] <-
      log(t(mapply(function(tauMin, tauRange) diff(seq(tauMin, tauMin+tauRange, length.out=nRatings-1)),
                   temp$tauMin, temp$tauRange)))

    inits[,(nCond + nRatings*2)] <- log(temp$m)

    logL <- apply(inits, MARGIN = 1,
                  function(p) try(ll_MratioF(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond), silent = TRUE))
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    inits <- inits[1:nInits,]

    noFitYet <- TRUE
    for (i in 1:nInits){

      m <- try(optim(par =  inits[i,],
                     fn = ll_MratioF, gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                     control = list(maxit = 10^4, reltol = 10^-4)))

      if (is.list(m)){
        for(j in 2:nRestart){
          try(m <- optim(par = m$par,
                         fn = ll_MratioF, gr = NULL,
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
        # c(as.vector(fit$par[nCond+nRatings-1] - rev(cumsum(c(exp(fit$par[(nCond+1):(nCond+nRatings-2)]))))),  as.vector(fit$par[nCond+nRatings-1]))
        exp(fit$par[nCond + nRatings*2]) * as.vector(fit$par[nCond+nRatings]) -
        rev( cumsum(c(exp(fit$par[(nCond+1):(nCond+nRatings-1)]))))

      res[,paste("theta_plus.",1:(nRatings-1), sep="")] <-
        # c(as.vector(fit$par[nCond+nRatings+1]), as.vector(fit$par[nCond+nRatings+1]) + as.vector(cumsum(c(exp(fit$par[(nCond+nRatings+2):(nCond + nRatings*2-1)])))))
        exp(fit$par[nCond + nRatings*2]) * as.vector(fit$par[nCond+nRatings]) +
        cumsum(c(exp(fit$par[(nCond+nRatings+1):(nCond + nRatings*2-1)])))

      res$m <- exp(fit$par[nCond + nRatings*2])

      res$negLogLik <- fit$value
      res$N <- N
      res$k <- k
      res$BIC <-  2 * fit$value + k * log(N)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }

