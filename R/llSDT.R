llSDT <-
function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p)) 
    locA <- -p[1:nCond]/2 
    locB <- p[1:nCond]/2 
    
    c_RA <- c(-Inf, p[nCond+nRatings] - 
                rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-1)])))),  
              p[nCond+nRatings])
    
    c_RB <- c(p[nCond+nRatings], p[nCond+nRatings] + 
                cumsum(c(exp(p[(nCond+nRatings+1):(length(p))]))), Inf)
    
    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    
    P_SBRB_SDT <- Vectorize(function(j,i) pnorm(q=c_RB[i+1], locB[j]) - pnorm(q=c_RB[i], locB[j]))
    P_SBRA_SDT <- Vectorize(function(j,i) pnorm(q=c_RA[i+1], locB[j]) - pnorm(q=c_RA[i], locB[j]))
    P_SARA_SDT <- Vectorize(function(j,i) pnorm(q=c_RA[i+1], locA[j]) - pnorm(q=c_RA[i], locA[j]))
    P_SARB_SDT <- Vectorize(function(j,i) pnorm(q=c_RB[i+1], locA[j]) - pnorm(q=c_RB[i], locA[j]))
    
    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SBRB_SDT) 
    p_SB_RA <- outer(1:nCond, 1:nRatings, P_SBRA_SDT)
    p_SA_RA <- outer(1:nCond, 1:nRatings, P_SARA_SDT)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SARB_SDT)
    
    p_SB_RB[(is.na(p_SB_RB))| is.nan(p_SB_RB)| p_SB_RB < 10^-64] <- 10^-64
    p_SB_RA[(is.na(p_SB_RA))| is.nan(p_SB_RA)| p_SB_RA < 10^-64] <- 10^-64
    p_SA_RB[(is.na(p_SA_RB))| is.nan(p_SA_RB)| p_SA_RB < 10^-64] <- 10^-64
    p_SA_RA[(is.na(p_SA_RA))| is.nan(p_SA_RA)| p_SA_RA < 10^-64] <- 10^-64
    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA, 
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))
    negLogL
  }
