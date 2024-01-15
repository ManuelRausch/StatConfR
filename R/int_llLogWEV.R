# in the logWEV model, the lognormal distribution should be shifted by theta otherwise the function we might get errors!
# also all rating criteria associated with R=1 should be greater than the disctimination criterion
# and all rating criteria associated with R=-1 should be smaller than the discrimination criterion!

ll_LogWEV <-
  function(p, N_SA_RA,N_SA_RB, N_SB_RA, N_SB_RB, nRatings, nCond){
    p <- c(t(p))
    ds <- cumsum(exp(p[1:(nCond)])) # enforce that sensitivity is ordered and all sensitivities are greater than 0
    locA <- -ds/2
    locB <- ds/2
    theta <- p[nCond+nRatings]

    c_RA <- c(-Inf, theta - # enforcde that all rating criteria associated with a are more negative than theta
                rev(cumsum(c(exp(p[(nCond+1):(nCond+nRatings-1)])))),
              theta)

    c_RB <- c(theta, theta + # enforcde that all rating criteria associated with B are more positive than theta
                cumsum(c(exp(p[(nCond+nRatings+1):(length(p))]))), Inf)

    sigma <- exp(p[nCond + nRatings*2])
    w <-  exp(p[nCond + nRatings*2+1])/(1+exp(p[nCond + nRatings*2+1]))

    p_SA_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SA_RB <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RA <- matrix(NA, nrow=nCond, ncol = nRatings)
    p_SB_RB <- matrix(NA, nrow=nCond, ncol = nRatings)

    P_SBRB <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (plnormALT(q = c_RB[i+1] - theta, (1 - w) * x + ds[j] * w - theta, sigma) - plnormALT(q = c_RB[i] - theta, (1 - w) * x + ds[j] * w - theta,  sigma)),
        lower = theta,
        upper = Inf,
        rel.tol = 10^-8)$value
    })
    P_SBRA <-  Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locB[j]) * (plnormALT(q = theta - c_RA[i+1], theta - (1 - w) * x + ds[j] * w, sigma) - plnormALT(q=theta-c_RA[i], theta - (1 - w) * x + ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARA <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (plnormALT(q=theta - c_RA[i+1], theta - (1 - w) * x + ds[j] * w, sigma) - plnormALT(q= theta - c_RA[i], theta - (1 - w) * x + ds[j] * w, sigma)),
        lower = -Inf,
        upper = theta,
        rel.tol = 10^-8)$value
    })
    P_SARB <- Vectorize(function(j,i){
      integrate(
        function(x) dnorm(x, locA[j]) * (plnormALT(q=c_RB[i+1] - theta, (1 - w) * x + ds[j] * w - theta, sigma) - plnormALT(q=c_RB[i] - theta, (1 - w) *x + ds[j] * w - theta, sigma)),
        lower = theta,
        upper= Inf,
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



plnormALT <-
  function (x, mean = exp(1/2), sd = sqrt(exp(1) - 1)*exp(1/2), log = FALSE)
  {
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean = as.vector(mean),
                             sd = as.vector(sd))
    for (i in c("x", "mean", "sd")) assign(i, arg.mat[, i])
    na.index <- is_na_matrix(arg.mat)
    if (all(na.index))
      y <- rep(NA, length(x))
    else {
      if (any(c(mean[!na.index], sd[!na.index]) < .Machine$double.eps))
        stop("All values of 'mean' and 'sd' must be positive.")
      cv <- sd/mean
      sdlog <- sqrt(log(1 + cv^2))
      meanlog <- log(mean) - (sdlog^2)/2
      y <- plnorm(x = x, meanlog = meanlog, sdlog = sdlog,
                  log = log)
    }
    if (!is.null(names.x))
      names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
  }


