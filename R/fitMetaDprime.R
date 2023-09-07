######################################
### Code to compute meta-d' in R #####
######################################

# Manuel Rausch

# computeMetaDprime computes meta-d' based on the matlab code provided by Brain Maniscalco http://www.columbia.edu/~bsm2105/type2sdt
# the method is slightly adapted to increase precision when optimazation routines of R are used
# first a coarse grid search. Afterwards, maximum likelihood estimation with several restarts of the fitting routine

# Arguments:
#   ratings: a factor with levels corresponding to the rating categories,
#            ordered from low to high. At least four different categories are required, five are recommended.
#   stimulus: a factor, levels corresponding to stimulus identities
#   correct: a vector with 1 indicating correct responses and 0 incorrect responses
#   Model: which version of meta-d' should be used
#         "Maniscalco" (default) uses the version of meta-d' by Maniscalco, B., & Lau, H. (2012). A signal detection theoretic method for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422-430.
#         "Fleming" uses the version of meta-d' by Fleming, S. M. (2017). HMeta-d: Hierarchical Bayesian estimation of metacognitive efficiency from confidence ratings. Neuroscience of Consciousness, 1, 1-14. https://doi.org/10.1093/nc/nix007
#   CI: calculate 95% confidence intervals
#   control: list of control parameter for the fitting routine

# Returns: a data.frame
# type 1 Hit rate: Observed
#   dprime:  estimated sensitivity from objective discrimination responses
#   cprime: estimatedstandardized response criterion
#   c raw response criterion
#   metaD estimated sensitivity from rating data
#   Ratio : metacognitive efficiency

computeMetaDprime   <- function(ratings, stimulus, correct,
                                ModelVersion = "Maniscalco",
                                CI = FALSE,
                                control = list(nRestart = 3,
                                               nStartValues = 5)) {
  if(!is.factor(ratings)) {
    stop ("ratings should be a factor with at least four levels ")
  }

  if(!is.factor(stimulus )|| length(levels(stimulus)) != 2) {
    stop("stimulus should be a factor with 2 levels")
  }
  if(!all(correct %in% c(0,1))) stop("correct should be 1 or 0")
  if(!any(correct == 0)) stop("There should be at least one erroneous response")
  if(!any(correct == 1)) stop("There should be at least one correct response")

  # prepare data

  nRatings <-  length(levels(ratings))
  nCriteria <- nRatings * 2 - 1
  abj_f <- 1 /(nRatings*2)
  abs_corrects <-  table(ratings[correct == 1], stimulus[correct == 1]) + abj_f
  abs_errors <- table(ratings[correct == 0], stimulus[correct == 0]) + abj_f

  abs_S1 <- c(rev(abs_errors[,2]),abs_corrects[,2])
  ratingHrs <- qnorm(1 - cumsum(abs_S1)/sum(abs_S1))
  abs_S2 <-  c(rev(abs_corrects[,1]), abs_errors[,1] )
  ratingFrs <-  qnorm(1 - cumsum(abs_S2)/sum(abs_S2))
  finits <- is.finite(ratingHrs) & is.finite(ratingFrs)
  ratingHrs <- as.vector(ratingHrs[finits])
  ratingFrs <- as.vector(ratingFrs[finits])

  nC_rS1 <- rev(as.vector(abs_corrects[,1]))
  nI_rS1 <- rev(as.vector(abs_errors[,2]))
  nC_rS2 <- as.vector(abs_corrects[,2])
  nI_rS2 <- as.vector(abs_errors[,1])

  # compute type 1 parameters based on formulae

  dprime <- ratingHrs[nRatings] - ratingFrs[nRatings]
  cs <- (-.5 * ( ratingHrs  + ratingFrs))
  cprime <- cs[nRatings]/dprime

  # use a coarse grid search to identify the most promising starting parameters
  temp <- expand.grid(d = seq(0, 5, length.out = 10),
                      tauMin =  seq(.1,2,length.out=10),  # position of the most conservative confidence criterion related to stimulus A
                      tauRange = seq(0.5,5,length.out=10))  # range of rating criteria stimulus B  #  position of the most liberal confidence criterion with respect to thet

  inits <- data.frame(matrix(data=NA, nrow= nrow(temp), ncol = 1 + (nRatings-1)*2))
  inits[,1] <- temp$d #  qnorm((temp$d + 10)/ 20)
  inits[,2:(nRatings-1)] <-
    log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                 temp$tauRange)))
  inits[,nRatings:(nRatings+1)] <- rep(log(temp$tauMin),2)
  inits[,(nRatings+2):(nRatings*2-1)] <-
    log(t(mapply(function(tauRange) rep(tauRange/(nRatings-1), nRatings-2),
                 temp$tauRange)))

  if(ModelVersion == "Maniscalco"){
    logL <- apply(inits, MARGIN = 1,
                  function(p) try(negLoglMetaD(p, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, cprime), silent = TRUE))
  }
  if (ModelVersion == "Fleming"){
    logL <- apply(inits, MARGIN = 1,
                  function(p) try(negLoglFleming(p, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, cs[nRatings]), silent = TRUE))
  }

  logL <- as.numeric(logL)
  inits <- inits[order(logL),]
  noFitYet <- TRUE

  if(ModelVersion == "Maniscalco"){
    for (j in 1:control$nStartValues){
      m <- try(optim(par = inits[j,], f = negLoglMetaD, gr = NULL,
                     nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                     nRatings = nRatings, cprime = cprime,
                     control = list(maxit = 10^6, reltol = 10^-8)), silent=T)
      for(i in 2:control$nRestart){
        try(m <- try(optim(par = m$par, f = negLoglMetaD, gr = NULL,
                           nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                           nRatings = nRatings, cprime = cprime,
                           control = list(maxit = 10^6, reltol = 10^-8)), silent=T),
            silent=T)
      }
      if ((class(m) != "try-error")){
        if (noFitYet) {
          fit <- m
          noFitYet <- FALSE
        } else {
          if (m$value < fit$value) fit <- m
        }
      }
    }
  }

  if(ModelVersion == "Fleming"){
    for (j in 1:control$nStartValues){
      m <- try(optim(par = inits[j,], f = negLoglFleming, gr = NULL,
                     nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                     nRatings = nRatings, type1_c = cs[nRatings],
                     control = list(maxit = 10^6, reltol = 10^-8)), silent=T)
      for(i in 2:control$nRestart){
        try(m <- try(optim(par = m$par, f = negLoglFleming, gr = NULL,
                           nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                           nRatings = nRatings, type1_c = cs[nRatings],
                           control = list(maxit = 10^6, reltol = 10^-8)), silent=T),
            silent=T)
      }
      if ((class(m) != "try-error")){
        if (noFitYet) {
          fit <- m
          noFitYet <- FALSE
        } else {
          if (m$value < fit$value) fit <- m
        }
      }
    }
  }


  result <- data.frame(
    dprime = dprime,
    cprime = cprime,
    c = dprime * cprime,
    metaD = NA, Ratio = NA, ModelVersion = ModelVersion)

  if(exists("fit")){
    if(class(fit)=="list"){
      result$metaD = fit$par[1] # pnorm(fit$par[1])*20 - 10
      result$Ratio =  result$metaD / dprime

      if (CI){
        if(ModelVersion == "Fleming"){
          Hess <-  try(optimHess(fit$par,  negLoglFleming,
                                 nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                                 nRatings = nRatings, type1_c = cs[nRatings]))
        }
        if(ModelVersion == "Maniscalco"){
          Hess <-  try(optimHess(fit$par,  negLoglMetaD,
                                 nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2,
                                 nRatings = nRatings, cprime = cprime))
        }


        if(!("try-error" %in% class(Hess))){
          SEs <-  sqrt(diag(solve(Hess)))
          result$metaD.SE <- SEs[1]
          result$CI95_lower <- result$metaD + qnorm(.025) * result$metaD.SE
          result$CI95_upper <- result$metaD + qnorm(.975) * result$metaD.SE
        }

      }
    }
  }
  result
}

negLoglMetaD <- function(parameters, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, cprime){
  metadprime <- parameters[1] # pnorm(parameters[1])*20 - 10
  S1mu <- -metadprime/2
  S2mu <- metadprime/2
  meta_c <- metadprime*cprime
  t2_rS1 <- c(-Inf, meta_c - rev(cumsum(exp(parameters[2:nRatings]))), meta_c)
  t2_rS2 <- c(meta_c, meta_c + cumsum(exp(parameters[(nRatings+1):length(parameters)])), Inf)

  prC_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)],S1mu) -
                pnorm(t2_rS1[1:nRatings],S1mu)) /pnorm(meta_c,S1mu)
  prI_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)], S2mu) -
                pnorm(t2_rS1[1:nRatings],S2mu) ) / pnorm(meta_c,S2mu)
  prC_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S2mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)], S2mu))) /
    (1 - pnorm(meta_c, S2mu))
  prI_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S1mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)],S1mu))) /
    (1 - pnorm(meta_c, S1mu))

  logL <-  -sum(nC_rS1*log(prC_rS1),nI_rS1*log(prI_rS1),
                nC_rS2*log(prC_rS2),nI_rS2*log(prI_rS2))
  prC_rS1[is.na( prC_rS1) | is.nan(prC_rS1) | prC_rS1 < 10^-10] <- 10^-10
  prI_rS1[is.na(prI_rS1) | is.nan(prI_rS1) | prI_rS1 < 10^-10] <- 10^-10
  prC_rS2[is.na( prC_rS2) | is.nan( prC_rS2) |  prC_rS2 < 10^-10] <- 10^-10
  prI_rS2[is.na(prI_rS2) | is.nan(prI_rS2) | prI_rS2 < 10^-10] <- 10^-10
  logL
}

negLoglFleming <- function(parameters, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, type1_c){
  metadprime <- parameters[1] # pnorm(parameters[1])*20 - 10
  S1mu <- -metadprime/2
  S2mu <- metadprime/2
  t2_rS1 <- c(-Inf, type1_c - rev(cumsum(exp(parameters[2:nRatings]))), type1_c)
  t2_rS2 <- c(type1_c, type1_c + cumsum(exp(parameters[(nRatings+1):length(parameters)])), Inf)

  prC_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)],S1mu) -
                pnorm(t2_rS1[1:nRatings],S1mu)) /pnorm(type1_c,S1mu)
  prI_rS1 <- (pnorm(t2_rS1[2:(nRatings+1)], S2mu) -
                pnorm(t2_rS1[1:nRatings],S2mu) ) / pnorm(type1_c,S2mu)
  prC_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S2mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)], S2mu))) /
    (1 - pnorm(type1_c, S2mu))
  prI_rS2 <- ((1- pnorm(t2_rS2[1:nRatings],S1mu)) -
                (1- pnorm(t2_rS2[2:(nRatings+1)],S1mu))) /
    (1 - pnorm(type1_c, S1mu))

  logL <-  -sum(nC_rS1*log(prC_rS1),nI_rS1*log(prI_rS1),
                nC_rS2*log(prC_rS2),nI_rS2*log(prI_rS2))
  prC_rS1[is.na( prC_rS1) | is.nan(prC_rS1) | prC_rS1 < 10^-10] <- 10^-10
  prI_rS1[is.na(prI_rS1) | is.nan(prI_rS1) | prI_rS1 < 10^-10] <- 10^-10
  prC_rS2[is.na( prC_rS2) | is.nan( prC_rS2) |  prC_rS2 < 10^-10] <- 10^-10
  prI_rS2[is.na(prI_rS2) | is.nan(prI_rS2) | prI_rS2 < 10^-10] <- 10^-10
  logL
}



