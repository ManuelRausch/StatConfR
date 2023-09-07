ll_Mratio <-
  function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p))
    ds <- cumsum(exp(p[1:(nCond)])) # enforce that sensitivity is ordered
    theta <- p[nCond+nRatings]
    locA1 <- -ds/2
    locB1 <- ds/2

    m_ratio <- exp(p[nCond + nRatings*2])
    metads <- m_ratio * ds
    locA2 <- -metads/2
    locB2 <- metads/2
    meta_c <- m_ratio * theta
    c_RA <- c(-Inf,
              meta_c - rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-1)])))),
              meta_c)
    c_RB <- c(meta_c, meta_c +
                cumsum(c(exp(p[(nCond+nRatings+1):(nCond + nRatings*2-1)]))),
              Inf)

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SBRB <- Vectorize(function(j,i){
      (1 - pnorm(theta, locB1[j])) * (pnorm(c_RB[i+1], locB2[j]) - pnorm(c_RB[i], locB2[j])) / (1 - pnorm(meta_c, locB2[j]))
    })
    P_SBRA <- Vectorize(function(j,i){
      pnorm(theta, locB1[j]) * (pnorm(c_RA[i+1], locB2[j]) - pnorm(c_RA[i], locB2[j])) / pnorm(meta_c, locB2[j])
    })
    P_SARA <-  Vectorize(function(j,i){
      pnorm(theta, locA1[j]) *
        (pnorm(c_RA[i+1], locA2[j]) - pnorm(c_RA[i], locA2[j])) /
        pnorm(meta_c, locA2[j])
    })
    P_SARB <-  Vectorize(function(j,i){
      (1 - pnorm(theta, locA1[j])) *
        (pnorm(c_RB[i+1], locA2[j]) - pnorm(c_RB[i], locA2[j])) /
        (1 - pnorm(meta_c, locA2[j]))
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

ll_MratioF <-
  function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p))
    ds <- cumsum(exp(p[1:(nCond)])) # enforce that sensitivity is ordered
    theta <- p[nCond+nRatings]
    locA1 <- -ds/2
    locB1 <- ds/2

    m_ratio <- exp(p[nCond + nRatings*2])
    metads <- m_ratio * ds
    locA2 <- -metads/2
    locB2 <- metads/2
    meta_c <-  theta # this is the version of the model used vy Fleming (2017)
    c_RA <- c(-Inf,
              meta_c - rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-1)])))),
              meta_c)
    c_RB <- c(meta_c, meta_c +
                cumsum(c(exp(p[(nCond+nRatings+1):(nCond + nRatings*2-1)]))),
              Inf)

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SBRB <- Vectorize(function(j,i){
      (1 - pnorm(theta, locB1[j])) * (pnorm(c_RB[i+1], locB2[j]) - pnorm(c_RB[i], locB2[j])) / (1 - pnorm(meta_c, locB2[j]))
    })
    P_SBRA <- Vectorize(function(j,i){
      pnorm(theta, locB1[j]) * (pnorm(c_RA[i+1], locB2[j]) - pnorm(c_RA[i], locB2[j])) / pnorm(meta_c, locB2[j])
    })
    P_SARA <-  Vectorize(function(j,i){
      pnorm(theta, locA1[j]) *
        (pnorm(c_RA[i+1], locA2[j]) - pnorm(c_RA[i], locA2[j])) /
        pnorm(meta_c, locA2[j])
    })
    P_SARB <-  Vectorize(function(j,i){
      (1 - pnorm(theta, locA1[j])) *
        (pnorm(c_RB[i+1], locA2[j]) - pnorm(c_RB[i], locA2[j])) /
        (1 - pnorm(meta_c, locA2[j]))
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
