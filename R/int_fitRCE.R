###   Functions to fit the response congruent evidence model
### Model based on Peters et al. (2017)
### the original model is extended to include a choice bias parameter so the
### the cognitive architecture underlying the decision is equivalent to SDT.

fitRCE <-
  function(N_SA_RA, N_SA_RB, N_SB_RA, N_SB_RB,
           nInits, nRestart, nRatings, nCond, nTrials){

    # coarse grid search to find promising initial values

    temp <- expand.grid(maxD =  seq(1, 5, 1),
                        theta = seq(-1/2,1/2, 1/2),
                        tauMin =  c(.1, .3, 1),  # position of the most conservative confidence criteria with respect to theta
                        tauRange = c(0.5, 1, 1.5, 2.5, 3.5))  #  position of the most liberal confidence criterion with respect to theta

    # number of parameters:
    # nCond sensitivity parameters
    # 1 type 1 bias parameter
    # 2 * (nRatings - 1) confidence criteria

    inits <- data.frame(matrix(data=NA, nrow= nrow(temp), ncol = nCond + nRatings*2 - 1))
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


    logL <- apply(inits, MARGIN = 1,
                  function(p) try(ll_llRCE(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond), silent = TRUE))
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    inits <- inits[1:nInits,]

    noFitYet <- TRUE
    for (i in 1:nInits){
      m <- try(optim(par =  inits[i,],
                     fn = ll_RCE, gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                     control = list(maxit = 10^4, reltol = 10^-4)))

      if (!inherits(m, "try-error")){
        for(j in 2:nRestart){
          try(m <- optim(par = m$par,
                         fn = ll_RCE, gr = NULL,
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

      res$negLogLik <- fit$value
      res$N <- nTrials
      res$k <- k
      res$BIC <-  2 * fit$value + k * log(nTrials)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(nTrials-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }

ll_RCE <-
  function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p))
    ds <- cumsum(exp(p[1:(nCond)])) # enforce that sensitivity is ordered

    b <- p[nCond+nRatings]
    theta <- p[nCond+nRatings]
    c_RA <- c(-Inf, p[nCond+nRatings-1] -
                rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-2)])))),
              p[nCond+nRatings-1], Inf)
    c_RB <- c(-Inf, p[nCond+nRatings+1], p[nCond+nRatings+1] +
                cumsum(c(exp(p[(nCond+nRatings+2):(nCond + nRatings*2-1)]))), Inf)
    sigma <- sqrt(1/2)

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SBRB <-  Vectorize(function(j, i){
      integrate(function(x) dnorm(x, mean=ds[j]/2, sigma) * pnorm(x-theta, mean=0, sigma),
                lower = c_RB[i],
                upper = c_RB[i+1],
                rel.tol = 10^-8)$value
    })
    P_SARB <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean = 0, sigma) * pnorm(x-theta, mean=ds[j]/2, sigma), # dnorm(x, mean=) * (1-pnorm(x, mean=b))
                lower = c_RB[i],
                upper = c_RB[i+1],
                rel.tol = 10^-8)$value
    })

    P_SBRA <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean= 0, sigma) * pnorm(x+theta, mean=ds[j]/2, sigma) , # dnorm(x, mean=ds[j]+b) * (1 - pnorm(x, mean=-b))
                lower = -c_RA[i+1], # braucht es hier einen Vorzeichenwechsel?
                upper = -c_RA[i],
                rel.tol = 10^-8)$value
    })
    P_SARA <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean = ds[j]/2, sigma) * pnorm(x+theta, mean=0, sigma),
                lower = -c_RA[i+1],
                upper = -c_RA[i],
                rel.tol = 10^-8)$value
    })

    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SBRB)
    p_SB_RA <- outer(1:nCond, 1:nRatings, P_SBRA)
    p_SA_RA <- outer(1:nCond, 1:nRatings, P_SARA)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SARB)

    p_SB_RB[(is.na(p_SB_RB))| is.nan(p_SB_RB)| p_SB_RB < 10^-64] <- 10^-64
    p_SB_RA[(is.na(p_SB_RA))| is.nan(p_SB_RA)| p_SB_RA < 10^-64] <- 10^-64
    p_SA_RB[(is.na(p_SA_RB))| is.nan(p_SA_RB)| p_SA_RB < 10^-64] <- 10^-64
    p_SA_RA[(is.na(p_SA_RA))| is.nan(p_SA_RA)| p_SA_RA < 10^-64] <- 10^-64

    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA,
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))

    negLogL
  }
