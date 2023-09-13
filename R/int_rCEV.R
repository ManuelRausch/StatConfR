rCEV <-
  function(p, n, nRatings, nConds){
    p <- c(t(p))
    locA <- -p[1:nConds]/2
    locB <- p[1:nConds]/2
    theta <- p[nConds+nRatings]
    c_RA <- c(-Inf, p[nConds+nRatings-1] -
                rev(cumsum(c(exp(p[(nConds+1):(nConds+nRatings-2)])))),
              p[nConds+nRatings-1], Inf)
    c_RB <- c(-Inf, p[nConds+nRatings+1], p[nConds+nRatings+1] +
                cumsum(c(exp(p[(nConds+nRatings+2):(nConds + nRatings*2-1)]))), Inf)
    sigma <- p[nConds + nRatings*2]
    w <-  p[nConds + nRatings*2 + 1]
    ds <- p[1:nConds] - mean(p[1:nConds])

    p_SA_RA <- matrix(NA, nrow=nConds, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nConds, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nConds, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nConds, ncol = nRatings)

    P_SBRB_CEV <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (pnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - pnorm(q=c_RB[i], (1 - w) *x + ds[j] * w, sigma)),
        lower = theta,
        upper = Inf,
        rel.tol = 10^-8)$value
    })
    P_SBRA_CEV <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (pnorm(q=c_RA[i+1], (1 - w) * x - ds[j] * w, sigma) - pnorm(q=c_RA[i], (1 - w) *x - ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARA_CEV <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (pnorm(q=c_RA[i+1], (1 - w) * x - ds[j] * w, sigma) - pnorm(q=c_RA[i], (1 - w) *x - ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARB_CEV <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (pnorm(q=c_RB[i+1], (1 - w) * x + ds[j] * w, sigma) - pnorm(q=c_RB[i], (1 - w) *x + ds[j] * w, sigma)),
        lower = theta,
        upper= Inf,
        rel.tol = 10^-8)$value
    })
    p_SB_RB <- outer(1:nConds, 1:nRatings, P_SBRB_CEV)
    p_SB_RA <- outer(1:nConds, 1:nRatings, P_SBRA_CEV)
    p_SA_RA <- outer(1:nConds, 1:nRatings, P_SARA_CEV)
    p_SA_RB <- outer(1:nConds, 1:nRatings, P_SARB_CEV)

    df <- expand.grid(n = 1:n,
                      condition = 1:nConds,
                      stimulus = c(-1, 1))
    response <- NULL
    probs <- rbind(cbind(p_SA_RA, p_SA_RB),
                   cbind(p_SB_RA, p_SB_RB))
    for (stimulus in 1:2) {
      for (condition in 1:nConds) {
        response <- c(response,
                      sample(x=1:(2*nRatings), size=n,
                             prob=probs[nConds*(stimulus-1)+condition,], replace = TRUE))
      }
    }
    df$response <- ifelse(response <= nRatings, -1, 1)
    df$correct <- as.numeric(df$stimulus==df$response)
    df$rating <- ifelse(df$response==1, response-nRatings, nRatings+1-response)
    return(df)
  }
