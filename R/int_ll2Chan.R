ll2Chan <-
function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p))
    ds <- exp(p[1:nCond])
    locA1 <- - ds /2
    locB1 <- ds /2
    metads <- ds * exp(p[nCond + nRatings*2])
    locA2 <- - metads/2
    locB2 <- metads/2
    theta <- p[nCond+nRatings]

    c_RA <- c(-Inf, p[nCond+nRatings-1] -
                rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-2)])))),
              p[nCond+nRatings-1], Inf)
    c_RB <- c(-Inf, p[nCond+nRatings+1], p[nCond+nRatings+1] +
                cumsum(c(exp(p[(nCond+nRatings+2):(nCond + nRatings*2-1)]))), Inf)

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SB_RB <- Vectorize(function(j,i){
      (1 - pnorm(theta, locB1[j])) * (pnorm(c_RB[i+1], locB2[j]) - pnorm(c_RB[i], locB2[j]))
    })
    P_SB_RA <- Vectorize(function(j,i){
      pnorm(theta, locB1[j]) * (pnorm(c_RA[i+1], locB2[j]) - pnorm(c_RA[i], locB2[j]))
    })
    P_SA_RA <-  Vectorize(function(j,i){
      pnorm(theta, locA1[j]) * (pnorm(c_RA[i+1], locA2[j]) - pnorm(c_RA[i], locA2[j]))
    })
    P_SA_RB <-  Vectorize(function(j,i){
      (1 - pnorm(theta, locA1[j])) * (pnorm(c_RB[i+1], locA2[j]) - pnorm(c_RB[i], locA2[j]))
    })

    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SB_RB)
    p_SB_RA <- outer(1:nCond, 1:nRatings, P_SB_RA)
    p_SA_RA <- outer(1:nCond, 1:nRatings, P_SA_RA)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SA_RB)

    p_SB_RB[(is.na(p_SB_RB))| is.nan(p_SB_RB)| p_SB_RB < 10^-64] <- 10^-64
    p_SB_RA[(is.na(p_SB_RA))| is.nan(p_SB_RA)| p_SB_RA < 10^-64] <- 10^-64
    p_SA_RB[(is.na(p_SA_RB))| is.nan(p_SA_RB)| p_SA_RB < 10^-64] <- 10^-64
    p_SA_RA[(is.na(p_SA_RA))| is.nan(p_SA_RA)| p_SA_RA < 10^-64] <- 10^-64

    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA,
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))
    negLogL
  }
