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
- metacognitive information thery
date: "31 October 2024"
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
- name: Sascha Meyen
  orcid: "0000-0001-6928-4126"
  affiliation: 3
- name: Sebastian Hellmann
  orcid: "0000-0001-6006-5103"
  equal-contrib: true
  affiliation: 2, 4
bibliography: paper.bib
affiliations:
- name: "Hochschule Rhein-Waal, Fakultät Gesellschaft und Ökonomie, Germany"
  index: 1
- name: "Katholische Universität Eichstätt-Ingolstadt, Philosophisch-pädagogische
    Fakultät, Germany"
  index: 2
- name: Eberhard Karls Universität Tübingen, Department of Computer Science, Tübingen, Germany
  index: 3
- name: TUM School of Management, Technische Universität München, Germany
  index: 4
---
  
# Summary
  
We present the `statConfR` package for R, which allows researchers to conveniently 
fit and compare nine different static models of decision confidence for 
binary discrimination tasks with confidence ratings: 
- the signal detection rating model [@Green1966], 
- the Gaussian noise model [@Maniscalco2016], 
- the independent Gaussian model [@Rausch2017], 
- the weighted evidence and visibility model [@Rausch2018], 
- the lognormal noise model [@Shekhar2020a],
- the lognormal weighted evidence and visibility model [@shekhar_how_2024], 
- the independent truncated Gaussian model [@rausch_measures_2023] based on the model specification 
used for the original meta-d$^\prime$/d$^\prime$ method [@Maniscalco2012; @Maniscalco2014], and 
- the independent truncated Gaussian model based on the model specification of Hmetad [@Fleming2017a]. 
In addition, the `statConfR` package provides functions for estimating different 
kinds of measures of metacognition: 
- meta-d$^\prime$/d$^\prime$, the most widely-used measure of metacognitive efficiency, allowing both @Maniscalco2012's and @Fleming2017a's model specification, 
- Information-theoretic measures [@dayan_metacognitive_2023], including
  - meta-I, an information-theoretic measures of metacognitive sensitivity, 
  - $meta-I_{1}^{r}$ and $meta-I_{2}^{r}$, two measures of metacognitive efficiency proposed by @dayan_metacognitive_2023, 
  - RMI, a novel measure of metacognitive accuracy, also derived from information theory.  
Finally, the `statConfR` package includes an example data set previously published in @hellmann_simultaneous_2023, with which the functions can be  tested. 
The `statConfR` reference manual provides documentation of each function of the latest release (https://cran.r-project.org/web/packages/statConfR/statConfR.pdf). 

# Statement of need

Cognitive models of confidence are currently used implicitly and explicitly in a wide range of research areas in the cognitive sciences: 
In perception research, confidence judgments can be used to quantify perceptual sensitivity based on receiver operating characteristics [@egan_operating_1959],
a method based on the signal detection rating model [@Green1966; @hautus_detection_2021]. 
In metacognition research, the most popular measure of metacognitive performance, 
the meta-d$^\prime$/d$^\prime$ method [@Maniscalco2012; @Maniscalco2014], 
implicitly relies on the independent truncated Gaussian model [@rausch_measures_2023]. 
Finally, confidence models have become a flourishing research topic in their own right 
[@boundy-singer_confidence_2022; @Desender2021; @guggenmos_reverse_2022; @hellmann_confidence_2024; @hellmann_simultaneous_2023; 
@pereira_evidence_2021; @Rausch2018; @Rausch2020; @Shekhar2020a; @shekhar_how_2024]. 
However, too few studies have empirically compared different confidence models [@Rausch2018; @Rausch2020; @rausch_measures_2023; 
@Shekhar2020a; @shekhar_how_2024], so there is still no consensus about the computational principles underlying confidence judgments 
[@rahnev_consensus_2022]. This is problematic because meta-d$^\prime$/d$^\prime$ can be biased by discrimination sensitivity, 
discrimination criteria, and/or confidence criteria if the generative model 
underlying the data is not the independent truncated Gaussian model [@rausch_measures_2023]. 
Likewise, receiver operating characteristics in rating experiments are only appropriate measures of discrimination sensitivity
if the assumptions of the signal detection rating model are correct [@Green1966; @hautus_detection_2021].
At the time of writing, `statConfR` is the only available package for an open software that allows
researchers to fit a set of static models of decision confidence. 
The ReMeta toolbox provides Python code to fit a variety of different confidence models [@guggenmos_reverse_2022], too, 
but some important models such as the independent truncated Gaussian model are missing. 
Previous studies modelling confidence have made their analysis scripts freely available on the OSF website 
[@Rausch2018; @Rausch2020; @rausch_measures_2023; @Shekhar2020a; @shekhar_how_2024], 
but these analysis scripts are often tailored to specific experiments and require time and effort to adapt to new experiments. 
In addition, the documentation of these scripts is not always sufficient to be used without
export knowledge in cognitive modelling. 
Finally, the lognormal noise model and the lognormal weighted evidence and visibility model were previously only available implemented in MATLAB, 
so `statConfR` makes these confidence models available to researchers who do not have access to MATLAB. 
The `statConfR` package also provides a faithful implementation of meta-d$^\prime$/d$^\prime$,
which has been originally implemented in MATLAB [@Maniscalco2012]. Fleming provides MATLAB and R code for Hmetad, 
a Bayesian hierarchical version of meta-d$^\prime$/d$^\prime$ [@Fleming2017a], 
but notably he specifies the model slightly differently as in the original meta-d$^\prime$/d$^\prime$ [@rausch_measures_2023]. 
To the best of our knowledge, there has been no open software available to estimate information-theoretic measures of metacognition up to now. 

An important limitation of the models implemented in `statConfR` is that the dynamics of the decision process are not taken into account.
This is a problem because confidence judgments are related to the dynamics of decision making [@hellmann_confidence_2024; @Pleskac2010; @Rahnev2020].
However, most previously proposed dynamical models of confidence do not include a parameter to represent metacognitive ability. 
There is one proposal for a dynamical measure of metacognitive efficiency, 
the v-ratio [@desender_dynamic_2022], which is based on two-stage signal detection theory [@Pleskac2010],
but two-stage signal detection theory has been outperformed by other models in a number of visual discrimination tasks 
[@hellmann_simultaneous_2023; @hellmann_confidence_2024; @shekhar_how_2024]. 
Thus, the static confidence models included in `statConfR` may still be useful for many researchers. 
 

## Contact

For comments, bug reports, and feature suggestions please contact 
<manuel.rausch@hochschule-rhein-waal.de> or <sebastian.hellmann@ku.de>
or [submit an issue](https://github.com/ManuelRausch/StatConfR/issues).

# Acknowledgements
    
This research was in part supported by grants RA2988/3-1, RA2988/4-1, and CRC 1233 (project 276693517) by the Deutsche Forschungsgemeinschaft. The funders had no role in software design, decision to publish, or preparation of the manuscript. Author contributions: Manuel Rausch: Conceptualization, Data curation, Funding acquisition, Software, Validation, Writing - original draft. Sascha Meyen: Software. Sebastian Hellmann: Conceptualization, Software, Writing - review and editing.
All participants who contributed to the data set provided written consent for participating in the experiment and for publishing their anonymized data in scientific repositories. The prodecure was approved by the Ethics Committee of the Katholische Universität Eichstätt-Ingolstadt. 

# References
