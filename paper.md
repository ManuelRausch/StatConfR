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
  
We present the statConfR package for R, which allows researchers to conveniently fit and compare nine different static models of decision confidence applicable to binary discrimination tasks with confidence ratings: the signal detection rating model [@Green1966], the Gaussian noise model[@Maniscalco2016], the independent Gaussian model [@Rausch2017], the weighted evidence and visibility model [@Rausch2018], the lognormal noise model [@Shekhar2020a], the lognormal weighted evidence and visibility model [@shekhar_how_2023], the independent truncated Gaussian model [@rausch_measures_2023] based on the meta-d′/d′ method [@Maniscalco2012; @Maniscalco2014], and the independent truncated Gaussian model based on the Hmetad method [@Fleming2017a]. In addition, the statConfR package provides functions to estimate meta-d′/d′, the most widely-used measure of metacognitive efficiency, allowing both @Maniscalco2012 and @Fleming2017a's model specification.

# Statement of need

Cognitive models of confidence are currently used implicitly and explicitly across a wide range of research areas in the Cognitive Sciences: In perception research, confidence judgements can be used to quantify perceptual sensitivity based on receiver operating characteristics [@egan_operating_1959], a method that relies on the signal detection rating model [@Green1966; @hautus_detection_2021]. In metacognition research, the most popular measure of metacognitive accuracy, the meta-d′/d′ method [@Maniscalco2012; @Maniscalco2014], depends on the independent truncated Gaussian model [@rausch_measures_2023]. Finally, models of confidence have become a flourishing research topic itself [@boundy-singer_confidence_2022; @Desender2021; @guggenmos_reverse_2022; @hellmann_confidence_2024; @hellmann_simultaneous_2023; @pereira_evidence_2021; @Rausch2018; @Rausch2020; @Shekhar2020a; @shekhar_how_2023]. However, up to date, too few studies have compared models of confidence empirically [@Rausch2018; @Rausch2020; @rausch_measures_2023; @Shekhar2020a, @shekhar_how_2023], which is why there is still no consensus about the computational principles underlying confidence judgments [@rahnev_consensus_2022]. This is problematic because meta-d′/d′ can be biased by discrimination sensitivity, discrimination criteria, and/or confidence criteria if the generative model underlying the data is not independent truncated Gaussian model [@rausch_measures_2023]. Likewise, receiver operating characteristics in rating experiments are appropriate measures of discrimination sensitivity only if the assumptions of the signal detection rating model are correct [@Green1966; @hautus_detection_2021].

At the time of writing this manuscript, statConfR is the only available package for an open software that allows researchers to fit a set of static models of decision confidence. The ReMeta-toolbox for MATLAB provides functions to fit a variety of different confidence models too [@guggenmos_reverse_2022], but several important models such as the independent truncated Gaussian model are missing. Previous studies modelling confidence have made their analysis scripts freely available at the osf website [REFERENCES MISSING], but these analysis scripts are often tailored to specific experiments and require time and effort to adapt to new experiments. In addition, the documentation of these scripts is not always sufficient to allow researchers without export knowledge in cognitive modelling to adapt these scripts for their own data sets. Finally, the lognormal noise model and the lognormal weighted evidence and visibility model have been previously  available only in MATLAB, which is why statConfR makes these confidence models available to researchers who do not have access to MATLAB.  
 
# Acknowledgements
    
This research was in part supported by grants RA2988/3-1 and RA2988/4-1 by the Deutsche Forschungsgemeinschaft. The funders had no role in study design, data collection, analysis, decision to publish, or preparation of the manuscript. We have no conflicts of interest to disclose.

# References
