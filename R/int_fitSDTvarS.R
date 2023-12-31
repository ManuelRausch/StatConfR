fitSDTvarS <-
function(ratings, stimulus, correct, condition){
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
    
    temp <- expand.grid(minD = c(0,.5),
                        maxD =  seq(3,14, length.out=4),
                        theta = seq(-1,1, 1),
                        tauMin =  seq(.1, 1, length.out=4),  # position of the most conservative confidence criteria with respect to theta
                        tauRange = seq(1,5, length.out=4),  #  position of the most liberal confidence criterion with respect to theta
                        a = c(1/4, 1/2,1, 2)) # range of stimulus and stimulus-unrelated variance at d=1
    
    inits <- data.frame(matrix(data=NA, nrow= nrow(temp), ncol = nCond + nRatings*2))
    inits[,1:nCond] <- t(mapply(function(minD, maxD) seq(minD, maxD, length.out = nCond), temp$minD, temp$maxD))
    inits[,(nCond+1):(nCond+nRatings-2)] <- 
      log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2), 
                   temp$tauRange)))
    inits[,nCond+(nRatings-1)] <- log(temp$tauMin)
    inits[,nCond+nRatings] <- temp$theta
    inits[,nCond+(nRatings+1)] <- log(temp$tauMin)
    inits[,(nCond+nRatings+2):(nCond + nRatings*2-1)] <- 
      log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2), 
                   temp$tauRange)))
    inits[,nCond + nRatings*2] <- log(temp$a)
    
    logL <- apply(inits, MARGIN = 1, 
                  function(p) try(llSDTvarS(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond),
                                  silent = TRUE))
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    inits <- inits[1:5,]
    
    nAttempts <- 25
    i <- 0
    noFitYet <- TRUE
    while (i < nAttempts |  noFitYet){
      i <- i+1
      m <- try(optim(par =  inits[i,], 
                     fn = llSDTvarS, gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB, 
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond, 
                     control = list(maxit = 10^6, reltol = 10^-8)))
      if (!inherits(m, "try-error")){
        inits <- rbind(inits, m$par)
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
      for (i in 1:nCond){
        res[paste("d",i, sep="")] <-  as.vector(fit$par[i])
      }
      res$theta <-  as.vector(fit$par[nCond+nRatings])
      res[,paste("cA",1:(nRatings-1), sep="")] <- 
        as.vector(fit$par[nCond+nRatings] - rev(cumsum(c(exp( fit$par[(nCond+1):(nCond+nRatings-1)])))))
      res[,paste("cB",1:(nRatings-1), sep="")] <- 
        as.vector(fit$par[nCond+nRatings] + cumsum(c(exp(fit$par[(nCond+nRatings+1):(nCond + nRatings*2-1)]))))
      
      res$a <- exp(fit$par[nCond + nRatings*2])
      
      res$negLogLik <- fit$value 
      res$N <- N
      res$k <- k 
      res$BIC <-  2 * fit$value + k * log(N)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }
