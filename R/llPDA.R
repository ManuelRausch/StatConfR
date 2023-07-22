llPDA <-
function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p)) 
    locA <- -p[1:nCond]/2 
    locB <- p[1:nCond]/2 
    theta <- p[nCond+nRatings]
    c_RA <- c(-Inf, p[nCond+nRatings-1] - 
                rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-2)])))),  
              p[nCond+nRatings-1], Inf)
    c_RB <- c(-Inf, p[nCond+nRatings+1], p[nCond+nRatings+1] + 
                cumsum(c(exp(p[(nCond+nRatings+2):(nCond + nRatings*2-1)]))), Inf)
    a <- exp(p[nCond + nRatings*2])
    sigma <- sqrt(a)
    
    p_SARA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SARB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SBRA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SBRB <- matrix(NA, nrow=nCond, ncol = nRatings)
    
    P_SBRB <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RB[i+1], x + locB[j]*a, sigma) - pnorm(q=c_RB[i], x + locB[j]*a, sigma)),
                lower = theta, 
                upper = Inf,
                rel.tol = 10^-8)$value
    })
    P_SBRA <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RA[i+1], x + locB[j]*a, sigma) - pnorm(q=c_RA[i], x + locB[j]*a, sigma)),
                lower = -Inf, 
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARA <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RA[i+1], x + locA[j]*a, sigma) - pnorm(q=c_RA[i], x + locA[j]*a, sigma)),
                lower = -Inf, 
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARB <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RB[i+1], x + locA[j]*a, sigma) - pnorm(q=c_RB[i], x +  locA[j]*a, sigma)),
                lower = theta, 
                upper= Inf,
                rel.tol = 10^-8)$value 
    })
    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SBRB) 
    p_SB_RA <- outer(1:nCond, 1:nRatings, P_SBRA)
    p_SA_RA <- outer(1:nCond, 1:nRatings, P_SARA)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SARB)
    
    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA, 
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))
    negLogL
  }
