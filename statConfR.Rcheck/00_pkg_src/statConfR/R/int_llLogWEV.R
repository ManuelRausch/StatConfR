ll_LogWEV <-
  function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){


    p <- c(t(p))
    ds <- cumsum(exp(p[1:(nCond)])) # enforce that sensitivity is ordered and all sensitivities are greater than 0
    locA <- -ds/2
    locB <- ds/2
    theta <- p[nCond+nRatings]

    c_RA <- c(0,
              -cumsum(c(exp(p[(nCond+1):(nCond+nRatings-1)]))),
              -Inf) # due to the lognormal distributions, the rating criteria are bounded by 0
    c_RB <- c(0, cumsum(c(exp(p[(nCond+nRatings+1):(nCond + nRatings*2-1)]))), Inf)

    sigma <- exp(p[nCond + nRatings*2])
    w <-  exp(p[nCond + nRatings*2+1])/(1+exp(p[nCond + nRatings*2+1]))

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SBRB <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (plnorm(q = c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - plnorm(q = c_RB[i], (1 - w) * x + ds[j] * w,  sigma)),
        lower = theta,
        upper = Inf,
        rel.tol = 10^-8)$value
    })
    P_SBRA <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (plnorm(q = -c_RA[i+1], -(1 - w) * x + ds[j] * w, sigma) - plnorm(q = -c_RA[i], -(1 - w) * x + ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARA <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (plnorm(q= -c_RA[i+1], -(1 - w) * x + ds[j] * w, sigma) - plnorm(q = -c_RA[i], -(1 - w) * x + ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARB <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (plnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - plnorm(q=c_RB[i], (1 - w) * x + ds[j] * w, sigma)),
        lower = theta,
        upper= Inf,
        rel.tol = 10^-8)$value
    })
    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SBRB)
    p_SB_RA <- outer(1:nCond, nRatings:1, P_SBRA)
    p_SA_RA <- outer(1:nCond, nRatings:1, P_SARA)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SARB)

    p_SB_RB[(is.na(p_SB_RB))| is.nan(p_SB_RB)| p_SB_RB < 10^-64] <- 10^-64
    p_SB_RA[(is.na(p_SB_RA))| is.nan(p_SB_RA)| p_SB_RA < 10^-64] <- 10^-64
    p_SA_RB[(is.na(p_SA_RB))| is.nan(p_SA_RB)| p_SA_RB < 10^-64] <- 10^-64
    p_SA_RA[(is.na(p_SA_RA))| is.nan(p_SA_RA)| p_SA_RA < 10^-64] <- 10^-64

    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA,
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))
    negLogL
  }

