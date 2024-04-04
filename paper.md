---
title: 'statConfR: An R Package for Static Models of Decision Confidence and Metacognition'
tags:
  - Cognitive modelling 
  - R code
  - signal detection theory
  - decision confidence
  - metacognition
  - meta-d′/d′
authors:
  - name: Manuel Rausch
    orcid: 0000-0002-5805-5544
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Sebastian Hellmann
    orcid: 0000-0001-6006-5103
    equal-contrib: true 
    affiliation: "2, 3"
affiliations:
 - name: Hochschule Rhein-Waal, Fakultät für Gesellschaft und Ökonomie, Germany
   index: 1
 - name: Katholische Universität Eichstätt-Ingolstadt, Philosophisch-pädagogische Fakultät, Germany
   index: 2
 - name: TUM School of Management, Technische Universität München, Germany
   index: 3
date: 04 April 2024
bibliography: paper.bib
---
  
# Summary
  
We present the statConfR package for R, which allows researchers to conveniently fit and compare nine different static models of decision confidence applicable to binary discrimination tasks with confidence ratings: the signal detection rating model [@Green1966], the Gaussian noise model[@Maniscalco2016], the independent Gaussian model [@Rausch2017], the weighted evidence and visibility model [@Rausch2018], the lognormal noise model [@Shekhar2020a], the lognormal weighted evidence and visibility model [@shekhar_how_2023], the independent truncated Gaussian model [@rausch_measures_2023] based on the meta-d′/d′ method [@Maniscalco2012; @Maniscalco2014], and the independent truncated Gaussian model based on the Hmetad method [@Fleming2017a]. In addition, the statConfR package provides functions for estimating meta-d′/d′, the most widely-used measure of metacognitive efficiency, allowing both @Maniscalco2012's and @Fleming2017a's model specification.

# Statement of need

Cognitive models of confidence are currently used implicitly and explicitly in a wide range of research areas in the cognitive sciences: In perception research, confidence judgments can be used to quantify perceptual sensitivity based on receiver operating characteristics [@egan_operating_1959], a method based on the signal detection rating model [@Green1966; @hautus_detection_2021]. In metacognition research, the most popular measure of metacognitive accuracy, the meta-d′/d′ method [@Maniscalco2012; @Maniscalco2014], implicitly relies on the independent truncated Gaussian model [@rausch_measures_2023]. Finally, confidence models have become a flourishing research topic in their own right [@boundy-singer_confidence_2022; @Desender2021; @guggenmos_reverse_2022; @hellmann_confidence_2024; @hellmann_simultaneous_2023; @pereira_evidence_2021; @Rausch2018; @Rausch2020; @Shekhar2020a; @shekhar_how_2023]. However, too few studies have empirically compared different confidence models [@Rausch2018; @Rausch2020; @rausch_measures_2023; @Shekhar2020a, @shekhar_how_2023], so there is still no consensus about the computational principles underlying confidence judgments [@rahnev_consensus_2022]. This is problematic because meta-d′/d′ can be biased by discrimination sensitivity, discrimination criteria, and/or confidence criteria if the generative model underlying the data is not the independent truncated Gaussian model [@rausch_measures_2023]. Likewise, receiver operating characteristics in rating experiments are only appropriate measures of discrimination sensitivity if the assumptions of the signal detection rating model are correct [@Green1966; @hautus_detection_2021].

At the time of writing, statConfR is the only available package for an open software that allows researchers to fit a set of static models of decision confidence. The ReMeta toolbox provides functions for MATLAB to also fit a variety of different confidence models [@guggenmos_reverse_2022], but some important models such as the independent truncated Gaussian model are missing. Previous studies modelling confidence have made their analysis scripts freely available on the OSF website [@rausch_full_2017; @rausch_full_2018; @rausch_full_2022; @shekhar_nature_2020-1; @shekhar_how_2022], but these analysis scripts are often tailored to specific experiments and require time and effort to adapt to new experiments. In addition, the documentation of these scripts is not always sufficient to be used without export knowledge in cognitive modelling. Finally, the lognormal noise model and the lognormal weighted evidence and visibility model were previously only available in MATLAB, so statConfR makes these confidence models available to researchers who do not have access to MATLAB.  
 
# Acknowledgements
    
This research was in part supported by grants RA2988/3-1 and RA2988/4-1 by the Deutsche Forschungsgemeinschaft. 

# References
