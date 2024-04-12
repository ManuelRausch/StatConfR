###   Functions to fit the WEV model
### Model based on Peters et al. (2017)
### the original model is extended to include a choice bias parameter so the
### the cognitive architecture underlying the decision is equivalent to SDT.

fitHeuris <-
  function(N_SA_RA, N_SA_RB, N_SB_RA, N_SB_RB,
           nInits, nRestart, nRatings, nCond){

     # coarse grid search to find promising initial values

    temp <- expand.grid(minD = seq(-.3, .8, length.out=4),
                        maxD =  c(.6, 1.2, 1.4, 1.9),
                        b = c(-.6, -.2, .2, .6),  # bias in favor of option "B
                        tauMin =  seq(-1.1, 2.4, length.out=4),  # position of the most conservative confidence criterion related to stimulus A
                        tauRange = seq(.1,3.8,  length.out=4)) # range of rating criteria stimulus B

    inits <- data.frame(matrix(data=NA, nrow= nrow(temp),
                               ncol = nCond + nRatings*2 -1))
    if(nCond==1)  {
      inits[,1] <-  log(temp$maxD)  }
    else{
      inits[,1:(nCond)] <-  log(t(mapply(function(maxD) diff(seq(0, maxD, length.out = nCond+1)), temp$maxD)))
    }
    inits[,nCond +1]  <- temp$tauMin
    inits[,(nCond+2):(nCond+nRatings-1)] <-
      log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                   temp$tauRange)))
    # inits[,nCond+(nRatings-1)] <- temp$tauMin
    inits[,nCond+nRatings] <- temp$b # theta
    inits[,nCond+(nRatings+1)] <- temp$tauMin
    inits[,(nCond+nRatings+2):(nCond + nRatings*2-1)] <-
      log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                   temp$tauRange)))

    logL <- apply(inits, MARGIN = 1,
                  function(p) try(llHeuris(p, N_SA_RA, N_SA_RB, N_SB_RA,N_SB_RB, nRatings, nCond), silent = TRUE))
    logL <- as.numeric(logL)
    nAttempts <- 5
    nRestart <- 4
    inits <- inits[order(logL),][1:nAttempts,]
    nIter <-  10^6
    noFitYet <- TRUE
    #print(paste("Initial grid search took...",as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2))," mins"))
    #print("Start fitting ... ")
    start <- inits[1,]

    try(fit <- optim(par =  start,
                     f = llHeuris(), gr = NULL,
                     N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                     N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                     control = list(maxit = 10^6, reltol = 10^-8)))
    if (!exists("fit") || class(fit) == "try-error"){
      i <- 2
      while(noFitYet && (i <= nAttempts)){
        start <- inits[i,]
        try(fit <- optim(par =  start,
                         f = llHeuris, gr = NULL,
                         N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                         N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                         control = list(maxit = 10^6, reltol = 10^-8)))
        if (exists("fit") && class(fit) == "list")  noFitYet <- FALSE
        i <- i+1
      }
    }

    if (exists("fit") && class(fit) == "list"){
      for (i in 1:nRestart){
        try(fit <- optim(par =  fit$par,
                         f = llHeuris, gr = NULL,
                         N_SA_RA = N_SA_RA,N_SA_RB = N_SA_RB,
                         N_SB_RA = N_SB_RA,N_SB_RB = N_SB_RB, nRatings = nRatings, nCond = nCond,
                         control = list(maxit = 10^6, reltol = 10^-8)))
      }
    }


    res <-  data.frame(matrix(nrow=1, ncol=0))
    if(class(fit) != "try-error"){
      k <- length(fit$par)
      N <- length(ratings)
      for (i in 1:nCond){
        res[paste("d",i, sep="")] <-  as.vector(fit$par[i])
      }
      res$b <-  as.vector(fit$par[nCond+nRatings])
      res[,paste("cA",1:(nRatings-1), sep="")] <-
        c(as.vector(fit$par[nCond+1]),
          as.vector(fit$par[nCond+1]) +
            as.vector(cumsum(c(exp(fit$par[(nCond+2):(nCond + nRatings-1)])))))

      res[,paste("cB",1:(nRatings-1), sep="")] <-
        c(as.vector(fit$par[nCond+nRatings+1]),
          as.vector(fit$par[nCond+nRatings+1]) +
            as.vector(cumsum(c(exp(fit$par[(nCond+nRatings+2):(nCond + nRatings*2-1)])))))


      res$negLogLik <- fit$value
      res$N <- N
      res$k <- k
      res$BIC <-  2 * fit$value + k * log(N)
      res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
      res$AIC <- 2 * fit$value + k * 2
    }
    res
  }

llHeuris <-
  function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p))

    ds <- p[1:nCond]
    b <- p[nCond+nRatings]
    c_RA <-  c(-Inf, p[nCond+1], p[nCond+1] +
                 cumsum(c(exp(p[(nCond+2):(nCond+nRatings-1)]))), Inf)
    # c_RA <- c(-Inf, p[nCond+nRatings-1],     rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-2)])))),       p[nCond+nRatings-1], Inf)
    c_RB <- c(-Inf, p[nCond+nRatings+1], p[nCond+nRatings+1] +
                cumsum(c(exp(p[(nCond+nRatings+2):(nCond + nRatings*2-1)]))), Inf)

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SBRB <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean=ds[j]+b) * pnorm(x, mean=-b),
                lower = c_RB[i],
                upper = c_RB[i+1],
                rel.tol = 10^-8)$value
    })
    P_SARB <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean = b) * pnorm(x, mean=ds[j]-b), # dnorm(x, mean=) * (1-pnorm(x, mean=b))
                lower = c_RB[i],
                upper = c_RB[i+1],
                rel.tol = 10^-8)$value
    })


    P_SBRA <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean=-b) * pnorm(x, mean=ds[j]+b) , # dnorm(x, mean=ds[j]+b) * (1 - pnorm(x, mean=-b))
                lower = c_RA[i],
                upper = c_RA[i+1],
                rel.tol = 10^-8)$value
    })
    P_SARA <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, mean= ds[j]-b) * pnorm(x, mean=b),
                lower = c_RA[i],
                upper = c_RA[i+1],
                rel.tol = 10^-8)$value
    })

    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SBRB)
    p_SB_RA <- outer(1:nCond, 1:nRatings, P_SBRA)
    p_SA_RA <- outer(1:nCond, 1:nRatings, P_SARA)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SARB)

    p_SB_RB[(is.na(p_SB_RB))| is.nan(p_SB_RB)| p_SB_RB < 1e-20] <- 1e-20
    p_SB_RA[(is.na(p_SB_RA))| is.nan(p_SB_RA)| p_SB_RA < 1e-20] <- 1e-20
    p_SA_RB[(is.na(p_SA_RB))| is.nan(p_SA_RB)| p_SA_RB < 1e-20] <- 1e-20
    p_SA_RA[(is.na(p_SA_RA))| is.nan(p_SA_RA)| p_SA_RA < 1e-20] <- 1e-20
    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA,
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))

    negLogL
  }
