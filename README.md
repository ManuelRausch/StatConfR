The `statConfR` package provides functions to fit static models of
decision-making and confidence derived from signal detection theory for
binary discrimination tasks, meta-d′/d′, the most prominent measure of
metacognitive efficiency, meta-I, an information-theoretic measures of
metacognitive sensitivity, as well as $`meta-I_{1}^{r}`$ and
$`meta-I_{2}^{r}`$, two information-theoretic measures of metacognitive
efficiency.

Fitting models of confidence can be used to test the assumptions
underlying meta-d′/d′. Several static models of decision-making and
confidence include a metacognition parameter that may serve as an
alternative when the assumptions of meta-d′/d′ assuming the
corresponding model provides a better fit to the data. The following
models are included:

- signal detection rating model (Green & Swets, 1966),
- Gaussian noise model (Maniscalco & Lau, 2016),
- weighted evidence and visibility model (Rausch et al., 2018),
- post-decisional accumulation model (Rausch et al., 2018),
- independent Gaussian model (Rausch & Zehetleitner, 2017),
- independent truncated Gaussian model (the model underlying the
  meta-d′/d′ method, see Rausch et al., 2023),
- lognormal noise model (Shekhar & Rahnev, 2021), and
- lognormal weighted evidence and visibility model (Shekhar & Rahnev,
  2023).

## Mathematical description of implemented models of confidence

The models included in the statConfR package are all based on signal
detection theory (Green & Swets, 1966). It is assumed that participants
select a binary discrimination response $`R`$ about a stimulus $`S`$.
Both $`S`$ and $`R`$ can be either -1 or 1. $`R`$ is considered correct
if $`S=R`$. In addition, we assume that in the experiment, there are
$`K`$ different levels of stimulus discriminability, i.e. a physical
variable that makes the discrimination task easier or harder. For each
level of discriminability, the function fits a different discrimination
sensitivity parameter $`d_k`$. If there is more than one sensitivity
parameter, we assume that the sensitivity parameters are ordered such as
$`0 < d_1 < d_2 < ... < d_K`$. The models assume that the stimulus
generates normally distributed sensory evidence $`x`$ with mean
$`S\times d_k/2`$ and variance of 1. The sensory evidence $`x`$ is
compared to a decision criterion $`c`$ to generate a discrimination
response $`R`$, which is 1, if $`x`$ exceeds $`c`$ and -1 else. To
generate confidence, it is assumed that the confidence variable $`y`$ is
compared to another set of criteria $`\theta_{R,i}, i=1,2,...,L-1`$,
depending on the discrimination response $`R`$ to produce a $`L`$-step
discrete confidence response. The different models vary in how $`y`$ is
generated (see below). The following parameters are shared between all
models:

- sensitivity parameters $`d_1, ..., d_K`$ ($`K`$: number of difficulty
  levels),
- decision criterion $`c`$,
- confidence criterion $`\theta_{-1,1}, ..., \theta_{-1,L-1},
  \theta_{1,1},  ,...,\theta_{1,L-1}`$ ($`L`$: number of confidence
  categories available for confidence ratings).

### Signal detection rating model (SDT)

According to SDT, the same sample of sensory evidence is used to
generate response and confidence, i.e., $`y=x`$. The confidence criteria
associated with $`R=-1`$ are more negative than the decision criterion
$`c`$, whereas the confidence criteria associated with $`R=1`$ are more
positive than $`c`$.

### Gaussian noise model (GN)

Conceptually, the Gaussian noise model reflects the idea that confidence
is informed by the same sensory evidence as the task decision, but
confidence is affected by additive Gaussian noise. According to GN,
$`y`$ is subject to additive noise and assumed to be normally
distributed around the decision evidence value $`x`$ with a standard
deviation $`\sigma`$, which is an additional free parameter.

### Weighted evidence and visibility model (WEV)

Conceptually, the WEV model reflects the idea that the observer combines
evidence about decision-relevant features of the stimulus with the
strength of evidence about choice-irrelevant features to generate
confidence. For this purpose, WEV assumes that $`y`$ is normally
distributed with a mean of $`(1-w)\times x+w \times d_k\times R`$ and
standard deviation $`\sigma`$. The standard deviation quantifies the
amount of unsystematic variability contributing to confidence judgments
but not to the discrimination judgments. The parameter $`w`$ represents
the weight that is put on the choice-irrelevant features in the
confidence judgment. The parameters $`w`$ and $`\sigma`$ are free
parameters in addition to the set of shared parameters.

### Post-decisional accumulation model (PDA)

PDA represents the idea of on-going information accumulation after the
discrimination choice. The parameter $`a`$ indicates the amount of
additional accumulation. The confidence variable is normally distributed
with mean $`x+S\times d_k\times a`$ and variance $`a`$. The parameter
$`a`$ is fitted in addition to the shared parameters.

### Independent Gaussian model (IG)

According to IG, the information used for confidence judgments is
generated independently from the sensory evidence used for the task
decision. For this purpose, it is assumed that $`y`$ is sampled
independently from $`x`$. The variable $`y`$ is normally distributed
with a mean of $`a\times d_k`$ and variance of 1. The additional
parameter $`m`$ represents the amount of information available for
confidence judgment relative to amount of evidence available for the
discrimination decision and can be smaller as well as greater than 1.

### Independent truncated Gaussian model: HMetad-Version (ITGc)

Conceptually, the two ITG models just as IG are based on the idea that
the information used for confidence judgments is generated independently
from the sensory evidence used for the task decision. However, in
contrast to IG, the two ITG models also reflect a form of confirmation
bias in so far as it is not possible to collect information that
contradicts the original decision. According to the version of ITG
consistent with the HMetad-method (Fleming, 2017), $`y`$ is sampled
independently from $`x`$ from a truncated Gaussian distribution with a
location parameter of $`S\times d_k \times m/2`$ and a scale parameter
of 1. The Gaussian distribution of $`y`$ is truncated in a way that it
is impossible to sample evidence that contradicts the original decision:
If $`R = -1`$, the distribution is truncated to the right of $`c`$. If
$`R = 1`$, the distribution is truncated to the left of $`c`$. The
additional parameter $`m`$ represents metacognitive efficiency, i.e.,
the amount of information available for confidence judgments relative to
amount of evidence available for discrimination decisions and can be
smaller as well as greater than 1.

### Independent truncated Gaussian model: Meta-d’-Version (ITGcm)

According to the version of the ITG consistent with the original meta-d’
method (Maniscalco & Lau, 2012, 2014), $`y`$ is sampled independently
from $`x`$ from a truncated Gaussian distribution with a location
parameter of $`S\times d_k \times m/2`$ and a scale parameter of 1. If
$`R = -1`$, the distribution is truncated to the right of $`m\times c`$.
If $`R = 1`$, the distribution is truncated to the left of
$`m\times c`$. The additional parameter $`m`$ represents metacognitive
efficiency, i.e., the amount of information available for confidence
judgments relative to amount of evidence available for the
discrimination decision and can be smaller as well as greater than 1.

### Logistic noise model (logN)

According to logN, the same sample of sensory evidence is used to
generate response and confidence, i.e., $`y=x`$ just as in SDT. However,
according to logN, the confidence criteria are not assumed to be
constant, but instead they are affected by noise drawn from a lognormal
distribution. In each trial, $`\theta_{-1,i}`$ is given by
$`c -  \epsilon_i`$. Likewise, $`\theta_{1,i}`$ is given by
$`c + \epsilon_i`$. The noise $`\epsilon_i`$ is drawn from a lognormal
distribution with the location parameter
$`\mu_{R,i} = \log(\left| \mu_{\theta_{R,i}} - c\right|)- 0.5 \times \sigma^{2}`$,
and scale parameter $`\sigma`$. $`\sigma`$ is a free parameter designed
to quantify metacognitive ability. It is assumed that the criterion
noise is perfectly correlated across confidence criteria, ensuring that
the confidence criteria are always perfectly ordered. Because
$`\theta_{-1,1}`$, …, $`\theta_{-1,L-1}`$, $`\theta_{1,1}`$, …,
$`\theta_{1,L-1}`$ change from trial to trial, they are not estimated as
free parameters. Instead, we estimate the means of the confidence
criteria, i.e., $`\mu_{\theta_{-1,1}}, ...,
\mu_{\theta_{-1,L-1}}, \mu_{\theta_{1,1}}, ...  \mu_{\theta_{1,L-1}}`$,
as free parameters.

### Logistic weighted evidence and visibility model (logWEV)

The logWEV model is a combination of logN and WEV proposed by .
Conceptually, logWEV assumes that the observer combines evidence about
decision-relevant features of the stimulus with the strength of evidence
about choice-irrelevant features. The model also assumes that noise
affecting the confidence decision variable is lognormal. According to
logWEV, the confidence decision variable is $`y`$ is equal to R × y’.
The variable y’ is sampled from a lognormal distribution with a location
parameter of $`(1-w)\times x\times R + w \times d_k`$ and a scale
parameter of $`\sigma`$. The parameter $`\sigma`$ quantifies the amount
of unsystematic variability contributing to confidence judgments but not
to the discrimination judgments. The parameter $`w`$ represents the
weight that is put on the choice-irrelevant features in the confidence
judgment. The parameters $`w`$ and $`\sigma`$ are free parameters.

## Measures of metacognition

### Meta-d’/d’

The conceptual idea of meta-d′ is to quantify metacognition in terms of
sensitivity in a hypothetical signal detection rating model describing
the primary task, under the assumption that participants had perfect
access to the sensory evidence and were perfectly consistent in placing
their confidence criteria (Maniscalco & Lau, 2012, 2014). Using a signal
detection model describing the primary task to quantify metacognition
allows a direct comparison between metacognitive accuracy and
discrimination performance because both are measured on the same scale.
Meta-d′ can be compared against the estimate of the distance between the
two stimulus distributions estimated from discrimination responses,
which is referred to as d′: If meta-d′ equals d′, it means that
metacognitive accuracy is exactly as good as expected from
discrimination performance. If meta-d′ is lower than d′, it means that
metacognitive accuracy is not optimal. It can be shown that the implicit
model of confidence underlying the meta-d’/d’ method is identical to
different versions of the independent truncated Gaussian model (Rausch
et al., 2023), depending on whether the original model specification by
Maniscalco and Lau (2012) or alternatively the specification by Fleming
(2017) is used. We strongly recommend to test whether the independent
truncated Gaussian models are adequate descriptions of the data before
quantifying metacognitive efficiency with meta-d′/d′.

### Information-theoretic measures of metacognition

Dayan (2023) proposed several measures of metacognition based on
quantities of information theory.

- Meta-I is a measure of metacognitive sensitivity defined as the mutual
  information between the confidence and accuracy and is calculated as
  the transmitted information minus the minimal information given the
  accuracy:

``` math
meta-I = I(Y; \hat{Y}, C) - I(Y; \hat{Y})
```
This is equivalent to Dayan’s formulation where meta-I is the
information that confidences transmit about the correctness of a
response:

``` math
meta-I = I(Y = \hat{Y}; C)
```
 - Meta-$`I_{1}^{r}`$ is meta-I normalized by the value of meta-I
expected assuming a signal detection model (Green & Swets, 1966) with
Gaussian noise, based on calculating the sensitivity index d’:

``` math
meta-I_{1}^{r} = meta-I / meta-I(d')
```
 - Meta-$`I_{2}^{r}`$ is meta-I normalized by its theoretical upper
bound, which is the information entropy of accuracy, $`H(Y = \hat{Y})`$:

``` math
meta-I_{2}^{r} = meta-I / H(Y = \hat{Y})
```

Notably, Dayan (2023) pointed out that a liberal or conservative use of
the confidence levels will affected the mutual information and thus all
information-theoretic measures of metacognition.

## Installation

The latest released version of the package is available on CRAN via

    install.packages("statConfR")

The easiest way to install the development version is using `devtools`
and install from GitHub:

    devtools::install_github("ManuelRausch/StatConfR")

<!-- without any dots, the code chunk will be shown, but not executed -->

## Usage

### Example data set

The package includes a demo data set from a masked orientation
discrimination task with confidence judgments (Hellmann et al., 2023,
Exp. 1).

``` r
library(statConfR)
data("MaskOri")
head(MaskOri)
```

    ##   participant stimulus response correct rating diffCond trialNo
    ## 1           1        0        0       1      0      8.3       1
    ## 2           1       90        0       0      4    133.3       2
    ## 3           1        0        0       1      0     33.3       3
    ## 4           1       90        0       0      0     16.7       4
    ## 5           1        0        0       1      3    133.3       5
    ## 6           1        0        0       1      0     16.7       6

### Fitting

The function `fitConfModels` allows the user to fit several confidence
models separately to the data of each participant. The data should be
provided via the argument `.data` in the form of a data.frame object
with the following variables in separate columns:

- stimulus (factor with 2 levels): The property of the stimulus which
  defines which response is correct
- diffCond (factor): The experimental manipulation that is expected to
  affect discrimination sensitivity
- correct (0-1): Indicating whether the choice was correct (1) or
  incorrect(0).
- rating (factor): A discrete variable encoding the decision confidence
  (high: very confident; low: less confident)
- participant (integer): giving the subject ID. The argument `model` is
  used to specify which model should be fitted, with ‘WEV’, ‘SDT’, ‘GN’,
  ‘PDA’, ‘IG’, ‘ITGc’, ‘ITGcm’, ‘logN’, and ‘logWEV’ as available
  options. If model=“all” (default), all implemented models will be fit,
  although this may take a while.

Setting the optional argument `.parallel=TRUE` parallizes model fitting
over all but 1 available core. Note that the fitting procedure takes may
take a considerable amount of time, especially when there are multiple
models, several difficulty conditions, and/or several confidence
categories. For example, if there are five difficulty conditions and
five confidence levels, fitting the WEV model to one single participant
may take 20-30 minutes on a 2.8GHz CPU. We recommend parallelization to
keep the required time tolerable.

``` r
fitted_pars <- fitConfModels(MaskOri, models=c("ITGcm", "WEV"), .parallel = TRUE) 
```

The output is then a data frame with one row for each combination of
participant and model and separate columns for each estimated parameter
as well as for different measures for goodness-of-fit (negative
log-likelihood, BIC, AIC and AICy<sub>c</sub>). These may be used for
statistical model comparisons.

``` r
head(fitted_pars)
```

    ##   model participant negLogLik   N  k      BIC     AICc      AIC          d_1
    ## 1 ITGcm           1 1046.6310 540 15 2187.636 2124.064 2123.262 5.454675e-10
    ## 2   WEV           1  984.7821 540 16 2070.229 2002.482 2001.564 7.878135e-02
    ## 3 ITGcm           2  887.1734 540 15 1868.720 1805.148 1804.347 3.829963e-02
    ## 4   WEV           2  854.9635 540 16 1810.592 1742.845 1741.927 1.694767e-01
    ## 5 ITGcm           3  729.1776 540 15 1552.729 1489.157 1488.355 2.268368e-01
    ## 6   WEV           3  720.8812 540 16 1542.428 1474.680 1473.762 4.514495e-01
    ##         d_2       d_3      d_4      d_5           c theta_minus.4 theta_minus.3
    ## 1 0.2202360 0.3548801 2.293571 3.425364 -0.08019056     -1.593619    -0.9176832
    ## 2 0.2627565 0.4221678 2.522239 2.990242 -0.10156005     -2.307768    -1.0284769
    ## 3 0.4414296 1.0067772 3.774056 4.826061 -0.32984521     -1.603116    -1.0359237
    ## 4 0.5410758 1.2167297 3.688688 4.533279 -0.38629823     -2.162288    -1.1900373
    ## 5 0.6150607 1.7070516 5.000353 6.558447 -0.44716359     -1.605777    -1.3117477
    ## 6 0.9393539 1.7224474 4.768763 6.022956 -0.48838134     -1.833939    -1.4803745
    ##   theta_minus.2 theta_minus.1 theta_plus.1 theta_plus.2 theta_plus.3
    ## 1    -0.5825042    -0.4093600    0.1129810   0.34925834    1.0131169
    ## 2    -0.2504396     0.2380873   -0.7998558  -0.06546803    1.3010215
    ## 3    -0.6044161    -0.4373389   -0.1826681   0.20172223    0.9977414
    ## 4    -0.1296185     0.6220593   -0.9230581   0.13155602    1.5182388
    ## 5    -0.8092707    -0.5866515   -0.2900164   0.05975546    0.9844860
    ## 6    -0.7322736    -0.1486786   -0.4187309   0.18464840    1.2401476
    ##   theta_plus.4        m     sigma         w         wAIC        wAICc
    ## 1     1.867748 1.143606        NA        NA 3.746813e-27 3.971060e-27
    ## 2     2.705993       NA 1.2098031 0.8788714 1.000000e+00 1.000000e+00
    ## 3     1.605164 1.036321        NA        NA 2.790743e-14 2.957770e-14
    ## 4     2.388335       NA 1.1305801 0.5395001 1.000000e+00 1.000000e+00
    ## 5     1.398138 1.076547        NA        NA 6.775184e-04 7.180389e-04
    ## 6     1.677248       NA 0.6994795 0.2729548 9.993225e-01 9.992820e-01
    ##           wBIC
    ## 1 3.203055e-26
    ## 2 1.000000e+00
    ## 3 2.385735e-13
    ## 4 1.000000e+00
    ## 5 5.762461e-03
    ## 6 9.942375e-01

It can be seen that the independent truncated Gaussian model is
consistently outperformed by the weighted evidence and visibility model,
which is why we would not recommend using meta-d′/d′ for this specific
task.

### Visualization

After obtaining model fits, it is strongly recommended to visualize the
prediction implied by the best fitting sets of parameters and to compare
the prediction with the actual data (Palminteri et al., 2017). The best
way to visualize the data is highly specific to the data set and
research question, which is why `statConfR` does not come with its own
visualization tools. This being said, here is an example for how a
visualization of model fit could look like:

<!-- Stuff where only the code should be shown and executed, but do not show R yapping  -->

``` r
library(tidyverse)
AggregatedData <- MaskOri %>%
  mutate(ratings = as.numeric(rating), diffCond = as.numeric(diffCond)) %>%
  group_by(participant, diffCond, correct ) %>% 
  dplyr::summarise(ratings=mean(ratings,na.rm=T)) %>%
  Rmisc::summarySEwithin(measurevar = "ratings",
                         withinvars = c("diffCond", "correct"), 
                         idvar = "participant",
                         na.rm = TRUE, .drop = TRUE) %>% 
  mutate(diffCond = as.numeric(diffCond))
AggregatedPrediction <- 
  rbind(fitted_pars %>%
          filter(model=="ITGcm") %>%
          group_by(participant) %>%
          simConf(model="ITGcm") %>% 
          mutate(model="ITGcm"), 
        fitted_pars %>%
          filter(model=="WEV") %>%
          group_by(participant) %>%
          simConf(model="WEV") %>% 
          mutate(model="WEV")) %>%
  mutate(ratings = as.numeric(rating) ) %>%
  group_by(participant, diffCond, correct, model ) %>% 
  dplyr::summarise(ratings=mean(ratings,na.rm=T)) %>%
  Rmisc::summarySEwithin(measurevar = "ratings",
                  withinvars = c("diffCond", "correct", "model"), 
                  idvar = "participant",
                  na.rm = TRUE, .drop = TRUE) %>% 
  mutate(diffCond = as.numeric(diffCond))
PlotMeans <- 
  ggplot(AggregatedPrediction, 
         aes(x = diffCond, y = ratings, color = correct)) + facet_grid(~ model) +
   ylim(c(1,5)) + 
   geom_line() +  ylab("confidence rating") + xlab("difficulty condition") +
   scale_color_manual(values = c("darkorange", "navy"),
                     labels = c("Error", "Correct response"), name = "model prediction") + 
  geom_errorbar(data = AggregatedData, 
                aes(ymin = ratings-se, ymax = ratings+se), color="black") + 
  geom_point(data = AggregatedData, aes(shape=correct), color="black") + 
  scale_shape_manual(values = c(15, 16),
                     labels = c("Error", "Correct response"), name = "observed data") + 
  theme_bw()
```

<!-- Show both the code and the output Figure!  -->

``` r
PlotMeans
```

<figure>
<img src="README_files/figure-gfm/unnamed-chunk-5-1.png"
alt="Predicted vs. observed confidence as a function of discriminability and correctness" />
<figcaption aria-hidden="true">Predicted vs. observed confidence as a
function of discriminability and correctness</figcaption>
</figure>

### Measuring metacognition

Assuming that the independent truncated Gaussian model provides a decent
account of the data (notably, this is not the case though in the demo
data set), the function `fitMetaDprime` can be used to estimate
meta-d′/d′ independently for each subject. The arguments `.data` and
`.parallel=TRUE` just in the same way the arguments of `fitConfModels`.
The argument `model` offers the user the choice between two model
specifications, either “ML” to use the original model specification used
by Maniscalco and Lau (2012, 2014) or “F” to use the model specification
by Fleming (2017)’s Hmetad method.

``` r
MetaDs <- fitMetaDprime(data = MaskOri, model="ML", .parallel = TRUE)
```

To estimate information theoretic measures of metacognition, we must
bring first the data into the correct format. Specifically, three
different types of inputs are accepted:

- A `data.frame` with variables “y” for true labels and “r” for
  confidence-binned responses. “y” needs to contain values -1 and +1
  while r needs to be a factor with ordered levels such that the first
  half of the levels correspond to predictions for y = -1 and the second
  half to predictions for y = +1.
- A counts `table` with joint absolute frequencies. Rows correspond to
  true labels (stimulus categories) and columns correspond to responses.
- A contingency `matrix` with joint relative frequencies (as before but
  normalized to sum up to 1).

``` r
OneSbj <- subset(MaskOri, participant == 1)
y <- ifelse(OneSbj$stimulus == 0, -1, 1)
r <- factor(ifelse(OneSbj$response == 0, -1, 1) * as.numeric(OneSbj$rating))
counts <- table(y, r)
```

Then, the different information-theoretic measures of metacognition can
be computed:

``` r
meta_I <- estimate_meta_I(counts)
meta_Ir1 <- estimate_meta_Ir1(counts)
meta_Ir1_acc <- estimate_meta_Ir1_acc(counts)
meta_Ir2 <- estimate_meta_Ir2(counts)
RMI <- estimate_RMI(counts)
```

### Documentation

The documentation of each function of the currently installed version of
`statConfR` can be accessed by typing *?functionname* into the console.

## Contributing to the package

The package is under active development. We are planning to implement
new models of decision confidence when they are published. Please feel
free to [contact us](malto::manuel.rausch@ku.de) to suggest new models
to implement in in the package, or to volunteer adding additional
models.

### Instruction for implementing custom models of decision confidence

**Only recommended for users with experience in cognitive modelling!**
For readers who want to use our open code to implement models of
confidence themselves, the following steps need to be taken:

- Derive the likelihood of a binary response ($`R=-1, 1`$) and a
  specific level of confidence ($`C=1,...K`$) according to the custom
  model and a set of parameters, given the binary stimulus
  ($`S=-1, 1`$), i.e. $`P(R, C | S)`$.
- Use one of the files named ‘int_ll*model*.R’ from the package sources
  and adapt the likelihood function according to your model. According
  to our convention, name the new file a ‘int_ll*yourmodelname*.R’. Note
  that all parameters are fitted on the reals, i.e. positive parameters
  should be transformed outside the log-likelihood function (e.g. using
  the logarithm) and back-transformed within the log-likelihood function
  (e.g. using the exponential).
- Use one of the files ‘int_fit*model*.R’ from the package sources and
  adapt the fitting function to reflect the new model.
  - The initial grid used during the grid search should include a
    plausible range of all parameters of your model.
  - If applicable, the parameters of the initial grid needs be
    transformed so the parameter vector for optimization is
    real-valued).
  - The optimization routine should call the new log-likelihood
    function.
  - If applicable, the parameter vector i obtained during optimization
    needs to back-transformation for the the output object `res`.
  - Name the new file according to the convention
    ‘int_fit*yourmodelname*.R’.
- Add your model and fitting-functions to the high-level functions
  `fitConf` and `fitConfModels`.
- Add a simulation function in the file ‘int_simulateConf.R’ which uses
  the same structure as the other functions but adapt the likelihood of
  the responses.

## Contact

For comments, bug reports, and feature suggestions please feel free to
write to either <manuel.rausch@hochschule-rhein-waal.de> or
<sebastian.hellmann@ku.de> or [submit an
issue](https://github.com/ManuelRausch/StatConfR/issues).

## References

- Dayan, P. (2023). Metacognitive Information Theory. Open Mind, 7,
  392–411. <https://doi.org/10.1162/opmi_a_00091>
- Fleming, S. M. (2017). HMeta-d: Hierarchical Bayesian estimation of
  metacognitive efficiency from confidence ratings. Neuroscience of
  Consciousness, 1, 1–14. <https://doi.org/10.1093/nc/nix007>
- Green, D. M., & Swets, J. A. (1966). Signal detection theory and
  psychophysics. Wiley.
- Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous
  modeling of choice, confidence, and response time in visual
  perception. Psychological Review, 130(6), 1521–1543.
  <https://doi.org/10.1037/rev0000411>
- Maniscalco, B., & Lau, H. (2012). A signal detection theoretic method
  for estimating metacognitive sensitivity from confidence ratings.
  Consciousness and Cognition, 21(1), 422–430.
  <https://doi.org/10.1016/j.concog.2011.09.021>
- Maniscalco, B., & Lau, H. (2016). The signal processing architecture
  underlying subjective reports of sensory awareness. Neuroscience of
  Consciousness, 1, 1–17. <https://doi.org/10.1093/nc/niw002>
- Maniscalco, B., & Lau, H. C. (2014). Signal Detection Theory Analysis
  of Type 1 and Type 2 Data: Meta-d, Response- Specific Meta-d, and the
  Unequal Variance SDT Model. In S. M. Fleming & C. D. Frith (Eds.), The
  Cognitive Neuroscience of Metacognition (pp. 25–66). Springer.
  <https://doi.org/10.1007/978-3-642-45190-4_3>
- Palminteri, S., Wyart, V., & Koechlin, E. (2017). The importance of
  falsification in computational cognitive modeling. Trends in Cognitive
  Sciences, 21(6), 425–433. <https://doi.org/10.1016/j.tics.2017.03.011>
- Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in
  masked orientation judgments is informed by both evidence and
  visibility. Attention, Perception, and Psychophysics, 80(1), 134–154.
  <https://doi.org/10.3758/s13414-017-1431-5>
- Rausch, M., & Zehetleitner, M. (2017). Should metacognition be
  measured by logistic regression? Consciousness and Cognition, 49,
  291–312. <https://doi.org/10.1016/j.concog.2017.02.007>
- Shekhar, M., & Rahnev, D. (2021). The Nature of Metacognitive
  Inefficiency in Perceptual Decision Making. Psychological Review,
  128(1), 45–70. <https://doi.org/10.1037/rev0000249>
- Shekhar, M., & Rahnev, D. (2024). How Do Humans Give Conﬁdence? A
  Comprehensive Comparison of Process Models of Perceptual
  Metacognition. Journal of Experimental Psychology: General, 153(3),
  656–688. <https://doi.org/10.1037/xge0001524>
