################################################################################
######### ANALYZE R-PACKAGE STATCONFR ##########################################
################################################################################

# Manuel Rausch, 15.04.2023

# To demonstrate that the package is working, we first fitted all models to the data of a
# masked orientation discrimination task, and then use the obtained parameter sets
# for a parameter recovery analysis. Unfortunately, running the complete script takes approximately
# on a PC with 15 cores. We upload the script anyway to document how we have varified the
# the functionality of our package.
# For those who want to test whether the package is working in principle, we
# provide a shorter script "QuickTests.R".
# We also uploaded the results of our own run of "TestScripts.R" to show the results of our
# parameter recovery analysis.

# 0) Preprare
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

# 3) Meta-d

# 1) Fit all models to the dataset from Hellmann et al. (2023) Exp. 1

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(statConfR)


head(MaskOri)

fitted_pars <-
  MaskOri %>%
  filter(participant < 3) %>% #uncomment this line if you don't have so much time for testing the package"
  fitConfModels( models = "all",
                 .parallel = TRUE, n.cores=9)


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
  geom_point(color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_SDT


# 2.2) GN

recov_pars_GN <-
  fitted_pars %>%
  filter(model=="GN") %>%
  group_by(participant) %>%
  simConf(model="GN") %>%
  fitConfModels(models = "GN", .parallel = TRUE)

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
  geom_point(color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_GN


# (iii) PDA

recov_pars_PDA <-
  fitted_pars %>%
  filter(model=="PDA") %>%
  group_by(participant) %>%
  simConf(model="PDA") %>%
  fitConfModels(models = "PDA", .parallel = TRUE)

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
  geom_point(color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_PDA


# (iv) IG

recov_pars_IG <-
  fitted_pars %>%
  filter(model=="IG") %>%
  group_by(participant) %>%
  simConf(model="IG") %>%
  fitConfModels(models = "IG", .parallel = TRUE)

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
  geom_point(color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_IG

# (v) WEV

recov_pars_WEV <-
  fitted_pars %>%
  filter(model=="WEV") %>%
  group_by(participant) %>%
  simConf(model="WEV") %>%
  fitConfModels(models = "WEV", .parallel = TRUE)

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
  geom_point(color="blue") +
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
  geom_point(color="blue") +
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
  geom_point(color="blue") +
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
  geom_point(color="blue") +
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
  geom_point(color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal()
Plot_recov_logWEV

# 3)

save(fitted_pars,
     recov_pars_SDT, Plot_recov_SDT,
     recov_pars_GN, Plot_recov_GN,
     recov_pars_logN, Plot_recov_logN,
     recov_pars_WEV, Plot_recov_WEV,
     recov_pars_logWEV, Plot_recov_logWEV,
     recov_pars_IG, Plot_recov_IG,
     recov_pars_ITGc, Plot_recov_ITGc,
     recov_pars_ITGcm, Plot_recov_ITGcm,
     recov_pars_PDA, Plot_recov_PDA,
     file = "TestResults.RData")
