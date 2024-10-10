# statConfR: Models of Decision Confidence and Measures of Metacognition

The statConfR package provides functions to fit static models of
decision-making and confidence derived from signal detection theory for
binary discrimination tasks, as well as meta-d′/d′ (Rausch & Hellmann,
2024). The package can be used to test the assumptions underlying
meta-d′/d′. Several models provide a metacognition parameter that may
serve as an alternative when the assumptions of meta-d′/d′ assuming the
corresponding model provides a better fit to the data. The following models are included:

- signal detection rating model,

- Gaussian noise model,

- weighted evidence and visibility model,

- post-decisional accumulation model,

- independent Gaussian model,

- independent truncated Gaussian model (the model underlying the
  meta-d′/d′ method, see Rausch et al., 2023),

- lognormal noise model, and

- lognormal weighted evidence and visibility model.

## Mathematical description of implemented models of confidence
The models included in the statConfR package are all based on signal detection theory (Green & Swets, 1966). 
It is assumed that participants select a binary discrimination response $R$ about a stimulus $S$.
Both $S$ and $R$ can be either -1 or 1. $R$ is considered correct if $S=R$.
In addition, we assume that in the experiment, there are $K$ different levels of stimulus discriminability, 
i.e. a physical variable that makes the discrimination task easier or harder.
For each level of discriminability, the function fits a different discrimination
sensitivity parameter $d_k$. If there is more than one sensitivity parameter,
we assume that the sensitivity parameters are ordered such as $0 < d_1 < d_2 < ... < d_K$.
The models assume that the stimulus generates normally distributed sensory evidence $x$ with mean $S\times d_k/2$
and variance of 1. The sensory evidence $x$ is compared to a decision criterion $c$ 
to generate a discrimination response $R$, which is 1, if $x$ exceeds $c$ and -1 else.
To generate confidence, it is assumed that the confidence variable $y$ is compared to another
set of criteria $\theta_{R,i}, i=1,2,...,L-1$, depending on the
discrimination response $R$ to produce a $L$-step discrete confidence response. The different models
vary in how $y$ is generated (see below). 
The parameters shared between all models are:
- sensitivity parameters $d_1, ..., d_K$ ($K$: number of difficulty levels),
- decision criterion $c$, 
- confidence criterion $\theta_{-1,1}, ..., \theta_{-1,L-1}, 
\theta_{1,1},  ,...,\theta_{1,L-1}$ ($L$: number of confidence categories available for confidence ratings).

### Signal detection rating model (SDT)
According to SDT, the same sample of sensory evidence is used to generate response and confidence, i.e.,
$y=x$. The confidence criteria associated with $R=-1$ are more negative than the 
decision criterion $c$, whereas the confidence criteria associated with $R=1$ are more positive than $c$ (Green & Swets, 1966).

### Gaussian noise model (GN)
According to GN, $y$ is subject to additive noise and assumed to be normally distributed
around the decision evidence value $x$ with a standard deviation $\sigma$, 
which is an additional free parameter (Maniscalco & Lau, 2016).

### Weighted evidence and visibility model (WEV)
WEV assumes that the observer combines evidence about decision-relevant features 
of the stimulus with the strength of evidence about choice-irrelevant features 
to generate confidence (Rausch et al., 2018). Thus, the WEV model assumes that $y$ is normally 
distributed with a mean of $(1-w)\times x+w \times d_k\times R$ and standard deviation $\sigma$.
The standard deviation quantifies the amount of unsystematic variability
contributing to confidence judgments but not to the discrimination judgments.
The parameter $w$ represents the weight that is put on the choice-irrelevant
features in the confidence judgment. The parameters $w$ and $\sigma$ are free parameters in
addition to the set of shared parameters.

### Post-decisional accumulation model (PDA)
PDA represents the idea of on-going information accumulation after the
discrimination choice (Rausch et al., 2018). The parameter $a$ indicates the amount of additional
accumulation. The confidence variable is normally distributed with mean
$x+S\times d_k\times a$ and variance $a$. The parameter $a$ is fitted in addition to the shared
parameters.

### Independent Gaussian model (IG)
According to IG, $y$ is sampled independently from $x$ (Rausch & Zehetleitner, 2017). 
The variable $y$ is normally distributed with a mean of $a\times d_k$ and variance
of 1. The additional parameter $m$ represents the amount of information available for confidence judgment
relative to amount of evidence available for the discrimination decision and can 
be smaller as well as greater than 1.

### Independent truncated Gaussian model: HMetad-Version (ITGc)
According to the version of ITG consistent with the HMetad-method 
(Fleming, 2017; see Rausch et al., 2023), $y$ is sampled independently
from $x$ from a truncated Gaussian distribution with a location parameter
of $S\times d_k \times m/2$ and a scale parameter of 1. The Gaussian distribution of $y$
is truncated in a way that it is impossible to sample evidence that contradicts
the original decision: If $R = -1$, the distribution is truncated to the
right of $c$. If $R = 1$, the distribution is truncated to the left
of $c$. The additional parameter $m$ represents metacognitive efficiency,
i.e., the amount of information available for confidence judgments relative to
amount of evidence available for discrimination decisions and  can be smaller
as well as greater than 1.

### Independent truncated Gaussian model: Meta-d'-Version (ITGcm)
According to the version of the ITG consistent with the original meta-d' 
method (Maniscalco & Lau, 2012, 2014; see Rausch et al., 2023),
$y$ is sampled independently from $x$ from a truncated Gaussian distribution with a location parameter
of $S\times d_k \times m/2$ and a scale parameter of 1. If $R = -1$, 
the distribution is truncated to the right of $m\times c$. If $R = 1$, 
the distribution is truncated to the left of  $m\times c$. The additional parameter $m$ 
represents metacognitive efficiency, i.e., the amount of information available for confidence judgments relative to
amount of evidence available for the discrimination decision and  can be smaller
as well as greater than 1.

### Logistic noise model (logN)
According to logN, the same sample of sensory evidence is used to generate response and confidence, i.e.,
$y=x$ just as in SDT (Shekhar & Rahnev, 2021). However, according to logN, the confidence criteria
are not assumed to be constant, but instead they are affected by noise drawn from
a lognormal distribution. In each trial, $\theta_{-1,i}$ is given
by $c -  \epsilon_i$. Likewise,  $\theta_{1,i}$ is given by
$c + \epsilon_i$. The noise $\epsilon_i$ is drawn from a lognormal distribution with
the location parameter $\mu_{R,i} = \log(\left| \theta_{R,i} - c\right|)- 0.5 \times \sigma^{2}$, 
$\mu_{R,i} = \log(\left| \bar{\theta}_{R,i} - c\right|)- 0.5 \times \sigma^{2}$. Continue here. 
$\mu_{R,i} =\log(\left|\overline{\theta}_{R,i}- c\right|) $ and
scale parameter $\sigma$. $\sigma$ is a free parameter designed to
quantify metacognitive ability. It is assumed that the criterion noise is perfectly
correlated across confidence criteria, ensuring that the confidence criteria
are always perfectly ordered. Because $\theta_{-1,1}$, ..., $\theta_{-1,L-1}$,
$\theta_{1,1}$, ..., $\theta_{1,L-1}$ change from trial to trial, they are not estimated
as free parameters. Instead, we estimate the means of the confidence criteria, i.e., $\overline{\theta}_{-1,1}, ...,
\overline{\theta}_{-1,L-1}, \overline{\theta}_{1,1}, ...  \overline{\theta}_{1,L-1}$,
as free parameters.

### \strong{Logistic Weighted Evidence and Visibility model (logWEV)
The logWEV model is a combination of logN and WEV, proposed by Shekhar and Rahnev (2023).
Conceptually, logWEV assumes that the observer combines evidence about decision-relevant features
of the stimulus with the strength of evidence about choice-irrelevant features (Rausch et al., 2018).
The model also assumes that noise affecting the confidence decision variable is lognormal
 in accordance with Shekhar and Rahnev (2021).
According to logWEV, the confidence decision variable is $y$ is equal to
$y^*\times R$. $y^*$ is sampled from a lognormal distribution with a location parameter
 of $(1-w)\times x\times R + w \times d_k$ and a scale parameter of $\sigma$.
 The parameter $\sigma$ quantifies the amount of unsystematic variability
contributing to confidence judgments but not to the discrimination judgments.
The parameter $w$ represents the weight that is put on the choice-irrelevant
features in the confidence judgment. $w$ and $\sigma$ are free parameters. 

## Installation

The latest released version of the package is available on CRAN via

    install.packages("statConfR")

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

## References

Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous
modeling of choice, confidence, and response time in visual perception.
Psychological Review. 130(6), 1521–1543.
[doi:10.1037/rev0000411](https://doi.org/10.1037/rev0000411)

Rausch, M., Hellmann, S. & Zehetleitner, M. (2023). Measures of
metacognitive efficiency across cognitive models of decision confidence.
Psychological Methods.
[doi:10.1037/met0000634](https://doi.org/10.1037/met0000634)

Rausch, M., & Hellmann, S. (2024). statConfR: An R Package for Static
Models of Decision Confidence and Metacognition. PsyArXiv.
[doi:10.31234/osf.io/dk6mr](https://doi.org/10.31234/osf.io/dk6mr)
