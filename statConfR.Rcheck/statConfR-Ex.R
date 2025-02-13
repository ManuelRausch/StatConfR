pkgname <- "statConfR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "statConfR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('statConfR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MaskOri")
### * MaskOri

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MaskOri
### Title: Data of 16 participants in a masked orientation discrimination
###   experiment (Hellmann et al., 2023, Exp. 1)
### Aliases: MaskOri
### Keywords: datasets

### ** Examples

data(MaskOri)
summary(MaskOri)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MaskOri", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimateMetaI")
### * estimateMetaI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimateMetaI
### Title: Estimate Measures of Metacognition from Information Theory
### Aliases: estimateMetaI

### ** Examples

# 1. Select two subjects from the masked orientation discrimination experiment
data <- subset(MaskOri, participant %in% c(1:2))
head(data)

# 2. Calculate meta-I measures with bias reduction (this may take 10 s per subject)
## No test: 
metaIMeasures <- estimateMetaI(data)
## End(No test)

# 3. Calculate meta-I measures for all participants without bias reduction (much faster)
metaIMeasures <- estimateMetaI(MaskOri, bias_reduction = FALSE)
metaIMeasures



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimateMetaI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitConf")
### * fitConf

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitConf
### Title: Fit a static confidence model to data
### Aliases: fitConf

### ** Examples

# 1. Select one subject from the masked orientation discrimination experiment
data <- subset(MaskOri, participant == 1)
head(data)

# 2. Use fitting function
## No test: 
  # Fitting takes some time (about 10 minutes on an 2.8GHz processor) to run:
  FitFirstSbjWEV <- fitConf(data, model="WEV")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitConf", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitConfModels")
### * fitConfModels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitConfModels
### Title: Fit several static confidence models to multiple participants
### Aliases: fitConfModels

### ** Examples

# 1. Select two subjects from the masked orientation discrimination experiment
data <- subset(MaskOri, participant %in% c(1:2))
head(data)

# 2. Fit some models to each subject of the masked orientation discrimination experiment
## No test: 
  # Fitting several models to several subjects takes quite some time
  # (about 10 minutes per model fit per participant on a 2.8GHz processor
  # with the default values of nInits and nRestart).
  # If you want to fit more than just two subjects,
  # we strongly recommend setting .parallel=TRUE
  Fits <- fitConfModels(data, models = c("SDT", "ITGc"), .parallel = FALSE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitConfModels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitMetaDprime")
### * fitMetaDprime

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitMetaDprime
### Title: title Compute measures of metacognitive sensitivity (meta-d')
###   and metacognitive efficiency(meta-d'/d') for data from one or several
###   subjects
### Aliases: fitMetaDprime

### ** Examples

# 1. Select two subject from the masked orientation discrimination experiment
data <- subset(MaskOri, participant %in% c(1:2))
head(data)

# 2. Fit meta-d/d for each subject in data
MetaDs <- fitMetaDprime(data, model="F", .parallel = FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitMetaDprime", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotConfModelFit")
### * plotConfModelFit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotConfModelFit
### Title: Plot the prediction of fitted parameters of one model of
###   confidence over the corresponding data
### Aliases: plotConfModelFit

### ** Examples

# 1. Fit some models to each subject of the masked orientation discrimination experiment
  # Normally, the fits should be created using the function fitConfModels
  # Fits <- fitConfModels(data, models = "WEV", .parallel = TRUE)
  # Here, we create the dataframe manually because fitting models takes about
  # 10 minutes per model fit per participant on a 2.8GHz processor.
  pars <- data.frame(participant = 1:16,
  d_1 = c(0.20, 0.05, 0.41, 0.03, 0.00, 0.01, 0.11, 0.03, 0.19, 0.08, 0.00,
  0.24, 0.00, 0.00, 0.25, 0.01),
  d_2 = c(0.61, 0.19, 0.86, 0.18, 0.17, 0.39, 0.69, 0.14, 0.45, 0.30, 0.00,
  0.27, 0.00, 0.05, 0.57, 0.23),
  d_3 = c(1.08, 1.04, 2.71, 2.27, 1.50, 1.21, 1.83, 0.80, 1.06, 0.68, 0.29,
  0.83, 0.77, 2.19, 1.93, 0.54),
  d_4 = c(3.47, 4.14, 6.92, 4.79, 3.72, 3.24, 4.55, 2.51, 3.78, 2.40, 1.95,
  2.55, 4.59, 4.27, 4.08, 1.80),
  d_5 = c(4.08, 5.29, 7.99, 5.31, 4.53, 4.66, 6.21, 4.67, 5.85, 3.39, 3.39,
  4.42, 6.48, 5.35, 5.28, 2.87),
  c = c(-0.30, -0.15, -1.37, 0.17, -0.12, -0.19, -0.12, 0.41, -0.27, 0.00,
  -0.19, -0.21, -0.91, -0.26, -0.20, 0.10),
  theta_minus.4 = c(-2.07, -2.04, -2.76, -2.32, -2.21, -2.33, -2.27, -2.29,
  -2.69, -3.80, -2.83, -1.74, -2.58, -3.09, -2.20, -1.57),
  theta_minus.3 = c(-1.25, -1.95, -1.92, -2.07, -1.62, -1.68, -2.04, -2.02,
  -1.84, -3.37, -1.89, -1.44, -2.31, -2.08, -1.53, -1.46),
  theta_minus.2 = c(-0.42, -1.40, -0.37, -1.96, -1.45, -1.27, -1.98, -1.66,
  -1.11, -2.69, -1.60, -1.25, -2.21, -1.68, -1.08, -1.17),
  theta_minus.1 = c(0.13, -0.90,  0.93, -1.71, -1.25, -0.59, -1.40, -1.00,
  -0.34, -1.65, -1.21, -0.76, -1.99, -0.92, -0.28, -0.99),
  theta_plus.1 = c(-0.62, 0.82, -2.77, 2.01, 1.39, 0.60, 1.51, 0.90, 0.18,
  1.62, 0.99,0.88, 1.67, 0.92, 0.18,  0.88),
  theta_plus.2 = c(0.15, 1.45, -1.13,2.17, 1.61, 1.24, 1.99, 1.55, 0.96, 2.44,
  1.53, 1.66, 2.00, 1.51, 1.08, 1.05),
  theta_plus.3 = c(1.40, 2.24, 0.77, 2.32, 1.80, 1.58, 2.19, 2.19, 1.54, 3.17,
  1.86, 1.85, 2.16, 2.09, 1.47, 1.70),
  theta_plus.4 = c(2.19, 2.40, 1.75, 2.58, 2.53, 2.24, 2.59, 2.55, 2.58, 3.85,
  2.87, 2.15, 2.51, 3.31, 2.27, 1.79),
  sigma = c(1.01, 0.64, 1.33, 0.39, 0.30, 0.75, 0.75, 1.07, 0.65, 0.29, 0.31,
  0.78, 0.39, 0.42, 0.69, 0.52),
  w = c(0.54, 0.50, 0.38, 0.38, 0.36, 0.44, 0.48, 0.48, 0.52, 0.46, 0.53, 0.48,
   0.29, 0.45, 0.51, 0.63))

# 2. Plot the predicted probabilities based on model and fitted parameters
  # against the observed relative frequencies.

  PlotFitWEV <- plotConfModelFit(MaskOri, pars, model="WEV")
  PlotFitWEV




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotConfModelFit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simConf")
### * simConf

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simConf
### Title: Simulate data according to a static model of confidence
### Aliases: simConf

### ** Examples

# 1. define some parameters
paramDf <- data.frame(d_1 = 0, d_2 = 2, d_3 = 4,c = .0,
theta_minus.2 = -2, theta_minus.1 = -1, theta_plus.1 = 1, theta_plus.2 = 2,
sigma = 1/2, w = 0.5, N = 500)
# 2. Simulate dataset
SimulatedData <- simConf(model = "WEV", paramDf)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simConf", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
