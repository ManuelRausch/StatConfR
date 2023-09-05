llCEVvarS <-
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
    sigma <- exp(p[nCond + nRatings*2])
    w <-  exp(p[nCond + nRatings*2 + 1])/(1+exp(p[nCond + nRatings*2 + 1]))
    ds <- p[1:nCond] - mean(p[1:nCond])
    sd <- #sqrt(exp(p[nCond+2])^2 + (p[1:nCond]/2)^2  )
      sqrt(1 + (p[1:nCond]/2)^2 * exp(p[nCond + nRatings*2 + 2])  )
    
    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    
    P_SBRB_CEV <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j], sd[j]) * (pnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - pnorm(q=c_RB[i], (1 - w) *x + ds[j] * w, sigma)),
        lower = theta, 
        upper = Inf,
        rel.tol = 10^-8)$value
    })
    P_SBRA_CEV <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j], sd[j]) * (pnorm(q=c_RA[i+1], (1 - w) * x - ds[j] * w, sigma) - pnorm(q=c_RA[i], (1 - w) *x - ds[j] * w, sigma)),
        lower = -Inf, 
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARA_CEV <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j], sd[j]) * (pnorm(q=c_RA[i+1], (1 - w) * x - ds[j] * w, sigma) - pnorm(q=c_RA[i], (1 - w) *x - ds[j] * w, sigma)),
        lower = -Inf, 
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARB_CEV <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j], sd[j]) * (pnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - pnorm(q=c_RB[i], (1 - w) *x + ds[j] * w, sigma)),
        lower = theta, 
        upper= Inf,
        rel.tol = 10^-8)$value 
    })
    p_SB_RB <- outer(1:nCond, 1:nRatings, P_SBRB_CEV) 
    p_SB_RA <- outer(1:nCond, 1:nRatings, P_SBRA_CEV)
    p_SA_RA <- outer(1:nCond, 1:nRatings, P_SARA_CEV)
    p_SA_RB <- outer(1:nCond, 1:nRatings, P_SARB_CEV)
    
    negLogL <- - sum (c(log(p_SB_RB) * N_SB_RB, log(p_SB_RA) * N_SB_RA, 
                        log(p_SA_RB) * N_SA_RB, log(p_SA_RA) * N_SA_RA))
    negLogL
  }
