################################################################################
######### ANALYZE R-PACKAGE STATCONFR ##########################################
################################################################################

# Manuel Rausch, 15.04.2023

# To demonstrate that the package is working, we first fitted all models to the data of a
# masked orientation discrimination task, and then use the obtained parameter sets
# for a parameter recovery analysis. Unfortunately, running the complete script takes several hours
# on a PC with 15 cores. We upload the script anyway to document how we have varified the
# the functionality of our package.
# For those who want to test whether the package is working in principle, we
# provide a shorter script "QuickTests.R".
# We also uploaded the results of our own run of "TestScripts.R" to show the results of our
# parameter recovery analysis.

# 0) Preparations
# 1) Fit all models to the dataset from Hellmann et al. (2023) Exp. 1
# 2) For each model: Parameter recovery based on simulated data using the fitted parameter sets
# 2.1) SDT
# 2.2) GN
# 2.3) PDA
# 2.4) IG
# 2.5) WEV
# 2.6) ITGc
# 2.7) ITGcm
# 2.8) logN
# 2.9) logWEV

# 3) meta-d′/d′
# 3.1) meta-d′/d′ using Maniscalco and Lau (2012)'s model specification
# 3.2) meta-d′/d′ using Fleming (2017)'s model specification

# 4. meta-I and co

# 0) Preparations

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # if you are not using R-Studio, specify here a working directory.

library(tidyverse)
library(statConfR)
if (file.exists("TestResults.RData")) load("TestResults.RData")

# 1) Fit all models to the dataset from Hellmann et al. (2023) Exp. 1


fitted_pars <-
  MaskOri %>%
  # filter(participant < 3) %>% #uncomment this line if you don't have so much time for testing the package"
  fitConfModels(models = "all",
                .parallel = TRUE)

PlotFitsBICWeights <-
  fitted_pars %>% #group_by(participant) %>%
  ggplot(aes(x=participant, y=wBIC, fill = model)) +
  geom_bar(stat="identity", color="black") +
  scale_x_continuous(breaks=unique(fitted_pars$participant)) +
  labs(fill = "Model")+
  ylab("Schwarz Weights") +
  theme_minimal()
PlotFitsBICWeights


# 2) Parameter recovery based on simulated data using the fitted parameter sets

# 2.1) SDT

recov_pars_SDT <-
  fitted_pars %>%
  filter(model=="SDT") %>%
  group_by(participant) %>%
  simConf(model="SDT") %>%
  fitConfModels(models = "SDT", .parallel = TRUE)

Plot_recov_SDT <-
  merge(recov_pars_SDT %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:theta_plus.4),
        fitted_pars %>%
          filter(model=="SDT") %>%
          select(participant, d_1:theta_plus.4) %>%
          pivot_longer(cols = d_1:theta_plus.4,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_SDT


# 2.2) GN

recov_pars_GN <-
  fitted_pars %>%
  filter(model=="GN") %>%
  group_by(participant) %>%
  simConf(model="GN") %>%
  fitConfModels(models = "GN",
                .parallel = TRUE)

Plot_recov_GN <-
  merge(recov_pars_GN %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:sigma),
        fitted_pars %>%
          filter(model=="GN") %>%
          select(participant, d_1:sigma) %>%
          pivot_longer(cols = d_1:sigma,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_GN


# (iii) PDA

recov_pars_PDA <-
  fitted_pars %>%
  filter(model=="PDA") %>%
  group_by(participant) %>%
  simConf(model="PDA") %>%
  fitConfModels(models = "PDA",
                .parallel = TRUE)

Plot_recov_PDA <-
  merge(recov_pars_PDA %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:b),
        fitted_pars %>%
          filter(model=="PDA") %>%
          select(participant, d_1:b) %>%
          pivot_longer(cols = d_1:b,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_PDA


# (iv) IG

recov_pars_IG <-
  fitted_pars %>%
  filter(model=="IG") %>%
  group_by(participant) %>%
  simConf(model="IG") %>%
  fitConfModels(models = "IG",
                .parallel = TRUE)

Plot_recov_IG <-
  merge(recov_pars_IG %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:m),
        fitted_pars %>%
          filter(model=="IG") %>%
          select(participant, d_1:m) %>%
          pivot_longer(cols = d_1:m,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_IG

# (v) WEV

recov_pars_WEV <-
  fitted_pars %>%
  filter(model=="WEV") %>%
  group_by(participant) %>%
  simConf(model="WEV") %>%
  fitConfModels(models = "WEV",
                .parallel = TRUE)

Plot_recov_WEV <-
  merge(recov_pars_WEV %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:w),
        fitted_pars %>%
          filter(model=="WEV") %>%
          select(participant, d_1:w) %>%
          pivot_longer(cols = d_1:w,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_WEV

# (vi) ITGc

recov_pars_ITGc <-
  fitted_pars %>%
  filter(model=="ITGc") %>%
  group_by(participant) %>%
  simConf(model="ITGc") %>%
  fitConfModels(models = "ITGc", .parallel = TRUE)

Plot_recov_ITGc <-
  merge(recov_pars_ITGc %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:m),
        fitted_pars %>%
          filter(model=="ITGc") %>%
          select(participant, d_1:m) %>%
          pivot_longer(cols = d_1:m,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_ITGc


# (vii) ITGcm

recov_pars_ITGcm <-
  fitted_pars %>%
  filter(model=="ITGcm") %>%
  group_by(participant) %>%
  simConf(model="ITGcm") %>%
  fitConfModels(models = "ITGcm", .parallel = TRUE)

Plot_recov_ITGcm <-
  merge(recov_pars_ITGc %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:m),
        fitted_pars %>%
          filter(model=="ITGcm") %>%
          select(participant, d_1:m) %>%
          pivot_longer(cols = d_1:m,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_ITGcm


# (viii) logN

recov_pars_logN <-
  fitted_pars %>%
  filter(model=="logN") %>%
  group_by(participant) %>%
  simConf(model="logN") %>%
  fitConfModels(models = "logN", .parallel = TRUE)

Plot_recov_logN <-
  merge(recov_pars_logN %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = c(d_1:c, M_theta_minus.4:M_theta_plus.4, sigma)),
        fitted_pars %>%
          filter(model=="logN") %>%
          select(participant, d_1:sigma) %>%
          pivot_longer(cols = d_1:sigma,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_logN


# (ix) logWEV

recov_pars_logWEV <-
  fitted_pars %>%
  filter(model=="logWEV") %>%
  group_by(participant) %>%
  simConf(model="logWEV") %>%
  fitConfModels( models = "logWEV", .parallel = TRUE)

Plot_recov_logWEV <-
  merge(recov_pars_logWEV %>%
          select(-model, -c(negLogLik:AIC)) %>%
          pivot_longer(cols = d_1:w),
        fitted_pars %>%
          filter(model=="logWEV") %>%
          select(participant, d_1:w) %>%
          pivot_longer(cols = d_1:w,
                       values_to = "true")) %>%
  ggplot(aes(x=true, y=value)) +
  facet_wrap(~ name, nrow=4, scales="free") + xlab("true parameter") + ylab("estimated parameter") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_logWEV

# 3) meta-d′/d′
# 3.1) meta-d′/d′ using Maniscalco and Lau (2012)'s model specification

recov_metaDprime_ML <-
  fitted_pars %>%
  filter(model=="ITGcm") %>%
  filter(participant!=11) %>%  # subject 11 performed very low.
  select(participant, d_3, c:theta_plus.4, m) %>%
  rename(d_1 = d_3) %>%
  mutate(N = 10000) %>% #
  group_by(participant) %>%
  simConf(model="ITGcm") %>%
  fitMetaDprime(model="ML", .parallel = TRUE)

Plot_recov_metaDprime_ML <-
  merge(recov_metaDprime_ML %>%
          select(participant, Ratio),
        fitted_pars %>%
          filter(model=="ITGcm") %>%
          select(participant, m)) %>%
  ggplot(aes(x=m, y=Ratio)) +  #scale_x_log10() + scale_y_log10() +
  xlab("m-parameter") + ylab("meta-d′/d′") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_metaDprime_ML


# 3.2) meta-d′/d′ using Fleming (2017)'s model specification

recov_metaDprime_F  <-
  fitted_pars %>%
  filter(model=="ITGc") %>%
  filter(participant!=11) %>%  # subject 11 performed very low, so meta-d'/d' will be unstable no matter whether the code works or not.
  select(participant, d_3, c:theta_plus.4, m) %>%
  rename(d_1 = d_3) %>%
  mutate(N = 10000) %>% # simulate 400 trials because 400 trials considered to be required to estimate meta-d′/d′
  group_by(participant) %>%
  simConf(model="ITGc") %>%
  fitMetaDprime(model="F", .parallel = TRUE)

Plot_recov_metaDprime_F <-
  merge(recov_metaDprime_F %>%
          select(participant, Ratio),
        fitted_pars %>%
          filter(model=="ITGc" ) %>%
          select(participant, m)) %>%
  ggplot(aes(x=m, y=Ratio)) + #scale_x_log10() + scale_y_log10() +
  xlab("m-parameter") + ylab("meta-d′/d′") +
  geom_point(color="purple") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_metaDprime_F


# 4. meta-I and co

MetaDs <- fitMetaDprime(data = MaskOri, model="ML", .parallel = TRUE)
MetaInfoMeasures  <- estimateMetaI(data = MaskOri, bias_reduction = F)


merge(MetaDs %>% select(participant, Ratio),
      MetaInfoMeasures %>%  select(participant, meta_Ir1)) %>%
  ggplot(aes(x=Ratio, y=meta_Ir1)) +
  geom_point() + geom_smooth(method="lm", se=F)+
  theme_minimal()


# 5) Plotting fits

PlotFitSDT <- plotConfModelFit(MaskOri, fitted_pars, model="SDT")
PlotFitGN <- plotConfModelFit(MaskOri, fitted_pars, model="GN")
PlotFitLogN <- plotConfModelFit(MaskOri, fitted_pars, model="logN")
PlotFitWEV <- plotConfModelFit(MaskOri, fitted_pars, model="WEV")
PlotFitLogWEV <- plotConfModelFit(MaskOri, fitted_pars, model="logWEV")
PlotFitITGcm <- plotConfModelFit(MaskOri, fitted_pars, model="ITGcm")
PlotFitITGc <- plotConfModelFit(MaskOri, fitted_pars, model="ITGc")
PlotFitIG <- plotConfModelFit(MaskOri, fitted_pars, model="IG")
PlotFitPDA <- plotConfModelFit(MaskOri, fitted_pars, model="PDA")


save(fitted_pars, PlotFitsBICWeights,
     recov_pars_SDT, Plot_recov_SDT,
     recov_pars_GN, Plot_recov_GN,
     recov_pars_logN, Plot_recov_logN,
     recov_pars_WEV, Plot_recov_WEV,
     recov_pars_logWEV, Plot_recov_logWEV,
     recov_pars_IG, Plot_recov_IG, # re-test that
     recov_pars_ITGc, Plot_recov_ITGc,
     recov_pars_ITGcm, Plot_recov_ITGcm,
     recov_pars_PDA, Plot_recov_PDA,

     recov_metaDprime_ML, Plot_recov_metaDprime_ML,
     recov_metaDprime_F, Plot_recov_metaDprime_F,
     MetaInfoMeasures,
     file = "TestResults.RData")
PlotFitSDT <- plotConfModelFit(MaskOri, fitatted_pars, model="SDT")
