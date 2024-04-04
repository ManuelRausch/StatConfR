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
  
We present the statConfR package for R, which allows researchers to conveniently fit and compare nine different static models of decision confidence applicable to binary discrimination tasks: the signal detection rating model [@green:1966], the Gaussian noise model [@Maniscalco:2016], the independent Gaussian model [@Rausch:2017], the weighted evidence and visibility model [@Rausch:2018], the lognormal noise model [@Shekhar:2021], the lognormal weighted evidence and visibility model [@Shekhar:2023], the independent truncated Gaussian model based on the original meta-d′/d′ method [@Maniscalco:2012; @Maniscalco:2014; @Rausch, 2023], and the independent truncated Gaussian model based on the Hmetad method  [@Fleming:2017; @Rausch, 2023]. In addition, the statConfR package provides functions to estimate meta-d′/d′, a widely-used measure of metacognitive efficiency based on both @Maniscalco:2012 and @Fleming:2017's model specification. 

# Statement of need

Cognitive models of confidence are currently used implicitly and explicitly across a wide range of research areas in the Cognitive Sciences: In perception research, confidence judgements can be used to quantify perceptual sensitivity based on receiver operating characteristics (Egan et al., 1959; Swets et al., 1961), a method that relies on the signal detection rating model [@Green:1966](Green & Swets, 1966; Hautus et al., 2021; Wickens, 2002). In metacognition research, the most popular measure of metacognitive accuracy, the meta-d′/d′ method (Maniscalco & Lau, 2012, 2014), depends on the independent truncated Gaussian model (Rausch et al., 2023). Finally, models of confidence have become a flourishing research topic itself (e.g. Boundy-Singer et al., 2022; Desender et al., 2021; Fleming & Daw, 2017; Guggenmos, 2022; Hellmann et al., 2023b, 2023a; Pereira et al., 2021; Rausch et al., 2018, 2020; Shekhar & Rahnev, 2021, 2023). However, the number of studies comparing model fit between different models of confidence is still relatively low (Rausch et al., 2018, 2020, 2023; Shekhar & Rahnev, 2021, 2023) and there is still no consensus about the computational principles underlying confidence judgments (Rahnev et al., 2022). This is problematic because meta-d′/d′ can be biased by discrimination sensitivity, discrimination criteria, and/or confidence criteria if the generative model underlying the data is not independent truncated Gaussian model (Rausch et al., 2023). Likewise, receiver operating characteristics in rating experiments are appropriate measures of discrimination sensitivity only if the assumptions of the signal detection rating model are correct (Green & Swets, 1966; Hautus et al., 2021). 

# Acknowledgements
    
This research was in part supported by grants RA2988/3-1 and RA2988/4-1 by the Deutsche Forschungsgemeinschaft. The funders had no role in study design, data collection, analysis, decision to publish, or preparation of the manuscript. We have no conflicts of interest to disclose.

# References
