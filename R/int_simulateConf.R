generateDataSDT <- function(paramDf){
  nCond <- sum(is.finite(c(t(paramDf[,grep(pattern = "d_", names(paramDf), value=T)]))))
  nRatings <- sum(is.finite(c(t(paramDf[,grep(pattern = "theta_minus", names(paramDf), value=T)]))))+1
  ds <- c(t(paramDf[,paste("d_", 1:nCond, sep="")]))
  locA <- -ds/2
  locB <- ds/2
  theta <- paramDf$c
  c_RA <- c(-Inf,c(t(paramDf[,paste("theta_minus.", (nRatings-1):1, sep="")])), theta)
  c_RB <- c(theta,c(t(paramDf[,paste("theta_plus.", 1:(nRatings-1), sep="")])), Inf)

  p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
  p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

  P_SBRB <- Vectorize(function(j,i) pnorm(q=c_RB[i+1], locB[j]) - pnorm(q=c_RB[i], locB[j]))
  P_SBRA <- Vectorize(function(j,i) pnorm(q=c_RA[i+1], locB[j]) - pnorm(q=c_RA[i], locB[j]))
  P_SARA <- Vectorize(function(j,i) pnorm(q=c_RA[i+1], locA[j]) - pnorm(q=c_RA[i], locA[j]))
  P_SARB <- Vectorize(function(j,i) pnorm(q=c_RB[i+1], locA[j]) - pnorm(q=c_RB[i], locA[j]))

  p_SB_RB <- cbind(stimulus = 1, response = 1,
                   mdply(.data=p_SB_RB, .fun= P_SBRB))
  p_SB_RA <- cbind(stimulus = 1, response = -1,
                   mdply(.data=p_SB_RA, .fun= P_SBRA))
  p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
  p_SA_RA <- cbind(stimulus = -1, response = -1,
                   mdply(.data=p_SA_RA, .fun= P_SARA))
  p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
  p_SA_RB <- cbind(stimulus = -1, response = 1,
                   mdply(.data=p_SA_RB, .fun= P_SARB))

  res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
  colnames(res) <- c("stimulus", "response", "condition",
                     "rating", "p")
  res$p[is.na(res$p) | is.nan(res$p)] <- 0
  nTrials <- round(paramDf$nTrials / 2 / nCond)

  f <- function(df){
    ind <- sample(x=1:(2*nRatings),size= nTrials,
                  prob = as.vector(df$p), replace=T)
    response <- df$response[ind]
    rating <- df$rating[ind]
    data.frame(response = response, rating = rating)
  }
  X <- ddply(res, .(stimulus, condition), f)
  X$condition <- factor(X$condition)
  X$correct <- 0
  X$correct[X$stimulus==X$response] <- 1
  X$rating <- factor(X$rating)
  X$stimulus <- factor(X$stimulus)
  X
}


generateDataNoisy <-
  function(paramDf){

    nCond <- length(grep(pattern = "d_", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "cA", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d", 1:nCond, sep="")]))
    locA <- -ds/2
    locB <- ds/2
    theta <- paramDf$theta
    c_RA <- c(-Inf,c(t(paramDf[,paste("cA", (nRatings-1):1, sep="")])), Inf)
    c_RB <- c(-Inf,c(t(paramDf[,paste("cB", 1:(nRatings-1), sep="")])), Inf)
    sigma <- paramDf$sigma

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

    P_SBRB_Noisy <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RB[i+1], x, sigma) - pnorm(q=c_RB[i], x, sigma)),
                lower = theta,
                upper = Inf,
                rel.tol = 10^-8)$value
    })
    P_SBRA_Noisy <-  Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locB[j]) * (pnorm(q=c_RA[i+1], x, sigma) - pnorm(q=c_RA[i], x, sigma)),
                lower = -Inf,
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARA_Noisy <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RA[i+1], x, sigma) - pnorm(q=c_RA[i], x, sigma)),
                lower = -Inf,
                upper = theta,
                rel.tol = 10^-8)$value
    })
    P_SARB_Noisy <- Vectorize(function(j,i){
      integrate(function(x) dnorm(x, locA[j]) * (pnorm(q=c_RB[i+1], x, sigma) - pnorm(q=c_RB[i], x, sigma)),
                lower = theta,
                upper= Inf,
                rel.tol = 10^-8)$value
    })

    p_SB_RB <- cbind(stimulus = "B", response = "B",
                     mdply(.data=p_SB_RB, .fun= P_SBRB_Noisy))
    p_SB_RA <- cbind(stimulus = "B", response = "A",
                     mdply(.data=p_SB_RA, .fun= P_SBRA_Noisy))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = "A", response = "A",
                     mdply(.data=p_SA_RA, .fun= P_SARA_Noisy))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = "A", response = "B",
                     mdply(.data=p_SA_RB, .fun= P_SARB_Noisy))

    res <- rbind(p_SB_RB,p_SB_RA,  p_SA_RA,p_SA_RB)
    colnames(res) <- c("stimulus", "response", "condition", "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    nTrials <- paramDf$nTrials

    f <- function(df){
      ind <- sample(x=1:(2*nRatings),size= nTrials, prob = as.vector(df$p), replace=T)
      response <- df$response[ind]
      rating <- df$rating[ind]
      data.frame(response = response, rating = rating)
    }
    X <- ddply(res, .(stimulus, condition), f)
    X$correct <- 0
    X$correct[X$stimulus==X$response] <- 1
    X$ratings <- factor(X$rating)
    X$stimulus <- factor(X$stimulus)
    X
  }

generateDataISDT <-
  function(paramDf){

    nCond <- length(grep(pattern = "d", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "cA", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d", 1:nCond, sep="")]))
    locA <- -ds/2
    locB <- ds/2
    theta <- paramDf$theta
    a <-  paramDf$a
    sigma <- sqrt(a)

    c_RA <- c(-Inf, c(t(paramDf[,paste("cA", (nRatings-1):1, sep="")])), Inf)
    c_RB <- c(-Inf, c(t(paramDf[,paste("cB", 1:(nRatings-1), sep="")])), Inf)

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

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

    p_SB_RB <- cbind(stimulus = "B", response = "B",
                     mdply(.data=p_SB_RB, .fun= P_SBRB))
    p_SB_RA <- cbind(stimulus = "B", response = "A",
                     mdply(.data=p_SB_RA, .fun= P_SBRA))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = "A", response = "A",
                     mdply(.data=p_SA_RA, .fun= P_SARA))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = "A", response = "B",
                     mdply(.data=p_SA_RB, .fun= P_SARB))

    res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
    colnames(res) <- c("stimulus", "response", "condition",
                       "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    nTrials <- paramDf$nTrials

    f <- function(df){
      ind <- sample(x=1:(2*nRatings),size= nTrials, prob = as.vector(df$p), replace=T)
      response <- df$response[ind]
      rating <- df$rating[ind]
      data.frame(response = response, rating = rating)
    }
    X <- ddply(res, .(stimulus, condition), f)
    X$correct <- 0
    X$correct[X$stimulus==X$response] <- 1
    X$ratings <- factor(X$rating)

    X$stimulus <- factor(X$stimulus)
    X
  }



generateDataWEV <-
  function(paramDf){
    require(plyr)
    nCond <- length(grep(pattern = "d", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "cA", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d", 1:nCond, sep="")]))
    locA <- -ds/2
    locB <- ds/2
    c_RA <- c(-Inf,c(t(paramDf[,paste("cA", (nRatings-1):1, sep="")])), Inf)
    c_RB <- c(-Inf,c(t(paramDf[,paste("cB", 1:(nRatings-1), sep="")])), Inf)

    theta <- paramDf$theta
    w <- paramDf$w
    sigma <- paramDf$sigma

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

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

    p_SB_RB <- cbind(stimulus = "B", response = "B",mdply(.data=p_SB_RB, .fun= P_SBRB_CEV))
    p_SB_RA <- cbind(stimulus = "B", response = "A",mdply(.data=p_SB_RA, .fun= P_SBRA_CEV))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = "A", response = "A",mdply(.data=p_SA_RA, .fun= P_SARA_CEV))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = "A", response = "B",mdply(.data=p_SA_RB, .fun= P_SARB_CEV))

    res <- rbind(p_SB_RB,p_SB_RA,  p_SA_RA,p_SA_RB)
    colnames(res) <- c("stimulus", "response", "condition", "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    nTrials <- paramDf$nTrials

    f <- function(df){
      ind <- sample(x=1:(2*nRatings),size= nTrials, prob = as.vector(df$p),
                    replace=T)
      response <- df$response[ind]
      rating <- df$rating[ind]
      data.frame(response = response, rating = rating)
    }
    X <- ddply(res, .(stimulus, condition), f)
    X$correct <- 0
    X$correct[X$stimulus==X$response] <- 1
    X$ratings <- factor(X$rating)
    X$stimulus <- factor(X$stimulus)
    X
  }



generateDataIndTruncML <-
  function(paramDf){

    nCond <- length(grep(pattern = "d", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "cA", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d", 1:nCond, sep="")]))
    locA1 <- -ds/2
    locB1 <- ds/2
    theta <- paramDf$theta

    m_ratio <- paramDf$m_ratio
    metads <- m_ratio * ds
    locA2 <- -metads/2
    locB2 <- metads/2
    meta_c <- m_ratio * theta

    c_RA <- c(-Inf, c(t(paramDf[,paste("cA", (nRatings-1):1, sep="")])), meta_c)
    c_RB <- c(meta_c, c(t(paramDf[,paste("cB", 1:(nRatings-1), sep="")])), Inf)

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

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

    p_SB_RB <- cbind(stimulus = "B", response = "B",
                     mdply(.data=p_SB_RB, .fun= P_SBRB))
    p_SB_RA <- cbind(stimulus = "B", response = "A",
                     mdply(.data=p_SB_RA, .fun= P_SBRA))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = "A", response = "A",
                     mdply(.data=p_SA_RA, .fun= P_SARA))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = "A", response = "B",
                     mdply(.data=p_SA_RB, .fun= P_SARB))

    res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
    colnames(res) <- c("stimulus", "response", "condition",
                       "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p) | res$p < 0] <- 0

    nTrials <- paramDf$nTrials

    f <- function(df){
      ind <- sample(x=1:(2*nRatings),size= nTrials, prob = as.vector(df$p), replace=T)
      response <- df$response[ind]
      rating <- df$rating[ind]
      data.frame(response = response, rating = rating)
    }
    X <- ddply(res, .(stimulus, condition), f)
    X$correct <- 0
    X$correct[X$stimulus==X$response] <- 1
    X$ratings <- factor(X$rating)
    X$stimulus <- factor(X$stimulus)
    X
  }

generateDataIndTruncF <-
  function(paramDf){

    nCond <- length(grep(pattern = "d", names(paramDf), value=T))
    nRatings <- (length(grep(pattern = "cA", names(paramDf), value=T)))+1
    ds <- c(t(paramDf[,paste("d", 1:nCond, sep="")]))
    locA1 <- -ds/2
    locB1 <- ds/2
    theta <- paramDf$theta

    m_ratio <- paramDf$m_ratio
    metads <- m_ratio * ds
    locA2 <- -metads/2
    locB2 <- metads/2
    meta_c <-  theta

    c_RA <- c(-Inf, c(t(paramDf[,paste("cA", (nRatings-1):1, sep="")])), meta_c)
    c_RB <- c(meta_c, c(t(paramDf[,paste("cB", 1:(nRatings-1), sep="")])), Inf)

    p_SA_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SA_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RA <- expand.grid(j = 1:nCond, i = 1:nRatings)
    p_SB_RB <- expand.grid(j = 1:nCond, i = 1:nRatings)

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

    p_SB_RB <- cbind(stimulus = "B", response = "B",
                     mdply(.data=p_SB_RB, .fun= P_SBRB))
    p_SB_RA <- cbind(stimulus = "B", response = "A",
                     mdply(.data=p_SB_RA, .fun= P_SBRA))
    p_SB_RA$i <- nRatings + 1 - p_SB_RA$i
    p_SA_RA <- cbind(stimulus = "A", response = "A",
                     mdply(.data=p_SA_RA, .fun= P_SARA))
    p_SA_RA$i <-  nRatings + 1 - p_SA_RA$i
    p_SA_RB <- cbind(stimulus = "A", response = "B",
                     mdply(.data=p_SA_RB, .fun= P_SARB))

    res <- rbind(p_SB_RB, p_SB_RA,  p_SA_RA, p_SA_RB)
    colnames(res) <- c("stimulus", "response", "condition",
                       "rating", "p")
    res$p[is.na(res$p) | is.nan(res$p)] <- 0
    nTrials <- paramDf$nTrials

    f <- function(df){
      ind <- sample(x=1:(2*nRatings),size= nTrials, prob = as.vector(df$p), replace=T)
      response <- df$response[ind]
      rating <- df$rating[ind]
      data.frame(response = response, rating = rating)
    }
    X <- ddply(res, .(stimulus, condition), f)
    X$correct <- 0
    X$correct[X$stimulus==X$response] <- 1
    X$ratings <- factor(X$rating)

    X$stimulus <- factor(X$stimulus)
    X
  }
