################################################################################
######### ANALYZE R-PACKAGE STATCONFR ##########################################
################################################################################

# Manuel Rausch, 15.04.2023
# We tested the full functionality of the package with the script "TestScript.R".
# However, running that script takes a considerable amount of time.
# QuickTests.R" is a script for those who want to test whether the package is working in principle, but
# do not have the time to run "TestScript.R".

# 1) Prepare

rm(list=ls())

library(tidyverse)


# Install and load the most recent development version of the package

if ("statConfR" %in% loadedNamespaces()){
  detach(package:statConfR, unload=TRUE)
}
devtools::install_github("ManuelRausch/StatConfR")

# Load the statConfR package

library(statConfR)

# load demo dataset

View(MaskOri)

# Fit a single confidence model to a single subject

FitOneSbj_SDT <-
  MaskOri %>%
  filter(participant == 1) %>%
  fitConf(model = "SDT")  # 20 s

# Fit multiple confidence models to a single subject with multiple cores
  # takes on my computer
  # may take twice as long if your computer has < 9 cores available.

system.time(FitOneSbj_allModels  <-
  MaskOri %>%
  filter(participant == 1) %>%
  fitConfModels(model = "all", .parallel=TRUE))


# Fit all nine confidence models to a single subject with multiple cores

FitSeveralSbj_MoreModels  <-
  MaskOri %>%
  filter(participant < 3) %>%  # In the script TestScript, we are doing this with all models and all subjects. However, that takes a while.
  fitConfModels(model = c("SDT", "ITGc"), .parallel=TRUE)   # 40 minutes


# Fit Meta-D prime
  # (not a good here because the ITGc is not a good fit to the data. This is for demonstration purposes only.)

MetaDs <- MaskOri %>%
  filter(diffCond == "33.3") %>%
  fitMetaDprime(model="ML", .parallel = TRUE)

