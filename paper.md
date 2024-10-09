---
title: 'statConfR: An R Package for Static Models of Decision Confidence and Measures
  of Metacognition'
tags:
- Cognitive modelling
- R code
- signal detection theory
- decision confidence
- metacognition
- "meta-d′/d′"
date: "08 October 2024"
output:
  pdf_document: default
  word_document: default
  html_document:
    df_print: paged
authors:
- name: Manuel Rausch
  orcid: "0000-0002-5805-5544"
  equal-contrib: true
  affiliation: 1, 2
- name: Sebastian Hellmann
  orcid: "0000-0001-6006-5103"
  equal-contrib: true
  affiliation: 2, 3
bibliography: paper.bib
affiliations:
- name: "Hochschule Rhein-Waal, Fakultät Gesellschaft und Ökonomie, Germany"
  index: 1
- name: "Katholische Universität Eichstätt-Ingolstadt, Philosophisch-pädagogische
    Fakultät, Germany"
  index: 2
- name: TUM School of Management, Technische Universität München, Germany
  index: 3
---
  
# Summary
  
We present the `statConfR` package for R, which allows researchers to conveniently fit and compare nine different static models of decision confidence applicable to binary discrimination tasks with confidence ratings: the signal detection rating model [@Green1966], the Gaussian noise model [@Maniscalco2016], the independent Gaussian model [@Rausch2017], the weighted evidence and visibility model [@Rausch2018], the lognormal noise model [@Shekhar2020a], the lognormal weighted evidence and visibility model [@shekhar_how_2024], the independent truncated Gaussian model [@rausch_measures_2023] based on the model specification used for the original meta-d$^\prime$/d$^\prime$ method [@Maniscalco2012; @Maniscalco2014], and the independent truncated Gaussian model based on the model specification of Hmetad [@Fleming2017a]. In addition, the `statConfR` package provides functions for estimating meta-d$^\prime$/d$^\prime$, the most widely-used measure of metacognitive efficiency, allowing both @Maniscalco2012's and @Fleming2017a's model specification. 
Finally, the `statConfR` package includes an example data set previously published in @hellmann_simultaneous_2023, with which the functions can be  tested. 

# Statement of need

Cognitive models of confidence are currently used implicitly and explicitly in a wide range of research areas in the cognitive sciences: In perception research, confidence judgments can be used to quantify perceptual sensitivity based on receiver operating characteristics [@egan_operating_1959], a method based on the signal detection rating model [@Green1966; @hautus_detection_2021]. In metacognition research, the most popular measure of metacognitive performance, the meta-d$^\prime$/d$^\prime$ method [@Maniscalco2012; @Maniscalco2014], implicitly relies on the independent truncated Gaussian model [@rausch_measures_2023]. Finally, confidence models have become a flourishing research topic in their own right [@boundy-singer_confidence_2022; @Desender2021; @guggenmos_reverse_2022; @hellmann_confidence_2024; @hellmann_simultaneous_2023; @pereira_evidence_2021; @Rausch2018; @Rausch2020; @Shekhar2020a; @shekhar_how_2024]. However, too few studies have empirically compared different confidence models [@Rausch2018; @Rausch2020; @rausch_measures_2023; @Shekhar2020a; @shekhar_how_2024], so there is still no consensus about the computational principles underlying confidence judgments [@rahnev_consensus_2022]. This is problematic because meta-d$^\prime$/d$^\prime$ can be biased by discrimination sensitivity, discrimination criteria, and/or confidence criteria if the generative model underlying the data is not the independent truncated Gaussian model [@rausch_measures_2023]. Likewise, receiver operating characteristics in rating experiments are only appropriate measures of discrimination sensitivity if the assumptions of the signal detection rating model are correct [@Green1966; @hautus_detection_2021].

At the time of writing, `statConfR` is the only available package for an open software that allows researchers to fit a set of static models of decision confidence. The ReMeta toolbox provides Python code to fit a variety of different confidence models [@guggenmos_reverse_2022], too, but some important models such as the independent truncated Gaussian model are missing. Previous studies modelling confidence have made their analysis scripts freely available on the OSF website [@Rausch2018; @Rausch2020; @rausch_measures_2023; @Shekhar2020a; @shekhar_how_2024], but these analysis scripts are often tailored to specific experiments and require time and effort to adapt to new experiments. In addition, the documentation of these scripts is not always sufficient to be used without export knowledge in cognitive modelling. Finally, the lognormal noise model and the lognormal weighted evidence and visibility model were previously only available implemented in MATLAB, so `statConfR` makes these confidence models available to researchers who do not have access to MATLAB. The `statConfR` package also provides a faithful implementation of meta-d$^\prime$/d$^\prime$, which has been originally implemented in MATLAB [@Maniscalco2012]. Fleming provides MATLAB and R code for Hmetad, a Bayesian hierarichical version of meta-d$^\prime$/d$^\prime$ [@Fleming2017a], but notably the model specification used for Hmetad is not the same as in meta-d$^\prime$/d$^\prime$ [@rausch_measures_2023].  

An important limitation of the models implemented in `statConfR` is that the dynamics of the decision process are not taken into account. This is a problem because confidence judgments are related to the dynamics of decision making [@hellmann_confidence_2024; @Pleskac2010; @Rahnev2020]. However, most previously proposed dynamical models of confidence do not include a parameter to represent metacognitive ability. There is one proposal for a dynamical measure of metacognitive efficiency, the v-ratio [@desender_dynamic_2022], which is based on two-stage signal detection theory [@Pleskac2010], but two-stage signal detection theory has been outperformed by other models in a number of visual discrimination tasks [@hellmann_simultaneous_2023; @hellmann_confidence_2024; @shekhar_how_2024]. Thus, the static confidence models included in `statConfR` may still be useful for many researchers. 
 
# Step-by-step example

## Installation

The latest released version of the package is available on CRAN via
```r
    install.packages("statConfR")
```
The easiest way to install the development version is using `devtools`
and install from GitHub:

    devtools::install_github("ManuelRausch/StatConfR")

## Usage

### Data structure

The package includes a demo data set from a masked orientation
discrimination task with confidence judgments (Hellmann et al., 2023,
Exp. 1.

``` r
library(statConfR)
data("MaskOri")
head(MaskOri)
```

    ##   participant stimulus correct rating diffCond trialNo
    ## 1           1        0       1      0      8.3       1
    ## 2           1       90       0      4    133.3       2
    ## 3           1        0       1      0     33.3       3
    ## 4           1       90       0      0     16.7       4
    ## 5           1        0       1      3    133.3       5
    ## 6           1        0       1      0     16.7       6

Data should be in the form of a data.frame object columns for following
variables:

- stimulus (factor with 2 levels): The property of the stimulus which
  defines which response is correct
- diffCond (factor): The experimental manipulation that is expected to
  affect discrimination sensitivity
- correct (0-1): Indicating whether the choice was correct (1) or
  incorrect(0).
- rating (factor): A discrete variable encoding the decision confidence
  (high: very confident; low: less confident)
- participant (integer): giving the subject ID.

### Fitting

It is strongly recommended that if metacognitive efficiency is to be
measured using the meta-d′/d′ method that researchers fist determine
whether the Independent Truncated Gaussian Model, the confidence model
implied by the meta-d′/d′ method, is an adequate description of the
data. Using the function fitConfModel, we can fit several confidence
models to the data of each participant. The argument
`.parallel=TRUE`allows for parallelization over all but one available
core.

    fitted_pars <- fitConfModels(MaskOri, models=c("SDT", "WEV"), .parallel = TRUE) 

This parallelizes the fitting process over participant-model
combinations. The output is then a data frame with one row for each
participant-model combination and columns for parameters and measures
for model performance (negative log-likelihood, BIC, AIC and AICc).
These may be used for quantitative model comparison.

``` r
head(fitted_pars)
```

    ##   model participant negLogLik    N  k      BIC     AICc      AIC    d_1    d_2
    ## 1   SDT           1  2721.256 1620 14 5545.975 5470.739 5470.513 0.0428 0.4593
    ## 2   WEV           1  2621.110 1620 16 5360.464 5274.520 5274.221 0.2027 0.6142
    ## 3   SDT           2  1946.258 1620 14 3995.979 3920.743 3920.517 0.0000 0.0950
    ## 4   WEV           2  1827.221 1620 16 3772.684 3686.741 3686.441 0.0512 0.1920
    ## 5   SDT           3  1706.178 1620 14 3515.818 3440.582 3440.356 0.2708 0.4673
    ## 6   WEV           3  1661.617 1620 16 3441.476 3355.533 3355.233 0.4146 0.8561
    ##      d_3    d_4    d_5       c theta_minus.4 theta_minus.3 theta_minus.2
    ## 1 1.0526 3.6806 4.7779 -0.2723       -1.5467       -1.0333       -0.6336
    ## 2 1.0797 3.4746 4.0799 -0.2957       -2.0665       -1.2485       -0.4152
    ## 3 0.8601 6.1410 8.0556 -0.1394       -2.0092       -1.9193       -1.4097
    ## 4 1.0412 4.1423 5.2886 -0.1475       -2.0441       -1.9500       -1.3982
    ## 5 1.9117 6.4257 7.5755 -1.1510       -1.9938       -1.6372       -1.2600
    ## 6 2.7115 6.9164 7.9863 -1.3743       -2.7625       -1.9192       -0.3724
    ##   theta_minus.1 theta_plus.1 theta_plus.2 theta_plus.3 theta_plus.4  sigma
    ## 1       -0.4543      -0.0944       0.2152       0.9850       1.5735     NA
    ## 2        0.1296      -0.6196       0.1544       1.3976       2.1879 1.0105
    ## 3       -0.9580       0.7857       1.3781       2.0879       2.2369     NA
    ## 4       -0.9030       0.8201       1.4484       2.2447       2.4030 0.6391
    ## 5       -1.1668      -1.1143      -0.7344       0.2961       0.9314     NA
    ## 6        0.9328      -2.7695      -1.1313       0.7714       1.7520 1.3289
    ##        w wAIC wAICc wBIC
    ## 1     NA    0     0    0
    ## 2 0.5361    1     1    1
    ## 3     NA    0     0    0
    ## 4 0.5020    1     1    1
    ## 5     NA    0     0    0
    ## 6 0.3818    1     1    1

If the Independent Truncated Gaussian model provides a decent account of
the data (which is not the case though in the demo dataset), it is
legitimate to quantify metacognitive efficiency with meta-d′/d′:

    MetaDs <- fitMetaDprime(subset(MaskOri, diffCond == "33.3"), 
                            model="ML", .parallel = TRUE)

## Contact

For comments, remarks, and questions please contact either
<manuel.rausch@ku.de> or <sebastian.hellmann@ku.de>
or [submit an issue](https://github.com/ManuelRausch/StatConfR/issues).

# Acknowledgements
    
This research was in part supported by grants RA2988/3-1 and RA2988/4-1 by the Deutsche Forschungsgemeinschaft. The funders had no role in software design, decision to publish, or preparation of the manuscript. Author contributions: Manuel Rausch: Conceptualization, Data curation, Funding acquisition, Software, Validation, Writing - original draft. Sebastian Hellmann: Conceptualization, Software, Writing - review and editing.
All participants who contributed to the dataset provided written consent for participating in the experiment and for publishing their anonymized data in scientific repositories. The prodecure was approved by the Ethics Committee of the Katholische Universität Eichstätt-Ingolstadt. 

# References
