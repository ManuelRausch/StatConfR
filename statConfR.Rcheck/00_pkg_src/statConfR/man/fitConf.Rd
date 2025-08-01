% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitConf.R
\name{fitConf}
\alias{fitConf}
\title{Fit a static confidence model to data}
\usage{
fitConf(data, model = "SDT", nInits = 5, nRestart = 4)
}
\arguments{
\item{data}{a \code{data.frame} where each row is one trial, containing following
variables:
\itemize{
\item \code{diffCond} (optional; different levels of discriminability,
should be a factor with levels ordered from hardest to easiest),
\item \code{rating} (discrete confidence judgments, should be a factor with levels ordered from lowest confidence to highest confidence;
otherwise will be transformed to factor with a warning),
\item \code{stimulus} (stimulus category in a binary choice task,
should be a factor with two levels, otherwise it will be transformed to
a factor with a warning),
\item \code{correct} (encoding whether the response was correct; should  be 0 for incorrect responses and 1 for correct responses)
}}

\item{model}{\code{character} of length 1. The generative model that should be
fitted. Models implemented so far: 'WEV', 'SDT', 'GN', 'PDA', 'IG',
'ITGc', 'ITGcm', 'logN', and 'logWEV'.}

\item{nInits}{\code{integer}. Number of starting values used for maximum likelihood optimization.
Defaults to 5.}

\item{nRestart}{\code{integer}. Number of times the optimization algorithm is restarted.
Defaults to 4.}
}
\value{
Gives data frame with one row and one column for each of the fitted parameters of the
selected model as well as additional information about the fit
(\code{negLogLik} (negative log-likelihood of the final set of parameters),
\code{k} (number of parameters), \code{N} (number of data rows),
\code{AIC} (Akaike Information Criterion; Akaike, 1974),
\code{BIC} (Bayes information criterion; Schwarz, 1978), and
\code{AICc} (AIC corrected for small samples; Burnham & Anderson, 2002))
}
\description{
The \code{fitConf} function fits the parameters of one static model of decision confidence,
provided by the \code{model} argument, to binary choices and confidence judgments.
See Details for the mathematical specification of the implemented models and
their parameters.
Parameters are fitted using a maximum likelihood estimation method with a
initial grid search to find promising starting values for the optimization.
In addition, several measures of model fit (negative log-likelihood, BIC, AIC, and AICc)
are computed, which can be used for a quantitative model evaluation.
}
\details{
The fitting routine first performs a coarse grid search to find promising
starting values for the maximum likelihood optimization procedure. Then the best \code{nInits}
parameter sets found by the grid search are used as the initial values for separate
runs of the Nelder-Mead algorithm implemented in \code{\link[stats]{optim}}.
Each run is restarted \code{nRestart} times.
\subsection{Mathematical description of models}{

The computational models are all based on signal detection theory (Green & Swets, 1966). It is assumed
that participants select a binary discrimination response \eqn{R} about a stimulus \eqn{S}.
Both \eqn{S} and \eqn{R} can be either -1 or 1.
\eqn{R} is considered correct if \eqn{S=R}.
In addition, we assume that there are \eqn{K} different levels of stimulus discriminability
in the experiment, i.e. a physical variable that makes the discrimination task easier or harder.
For each level of discriminability, the function fits a different discrimination
sensitivity parameter \eqn{d_k}. If there is more than one sensitivity parameter,
we assume that the sensitivity parameters are ordered such as \eqn{0 < d_1 <  ... < d_K}.
The models assume that the stimulus generates normally distributed sensory evidence \eqn{x} with mean \eqn{S\times d_k/2}
and variance of 1. The sensory evidence \eqn{x} is compared to a decision
criterion \eqn{c} to generate a discrimination response
\eqn{R}, which is 1, if \eqn{x} exceeds \eqn{c} and -1 else.
To generate confidence, it is assumed that the confidence variable \eqn{y} is compared to another
set of criteria \eqn{\theta_{R,i}, i = 1, ..., L-1}, depending on the
discrimination response \eqn{R} to produce a \eqn{L}-step discrete confidence response.
The number of thresholds will be inferred from the number of steps in the
\code{rating} column of \code{data}. Thus, the parameters shared between all models are:
\itemize{
\item sensitivity parameters \eqn{d_1},...,\eqn{d_K} (\eqn{K}: number of difficulty levels)
\item decision criterion \eqn{c}
\item confidence criterion \eqn{\theta_{-1,1}},\eqn{\theta_{-1,2}},
..., \eqn{\theta_{-1,L-1}}, \eqn{\theta_{1,1}},  \eqn{\theta_{1,2}},...,
\eqn{\theta_{1,L-1}} (\eqn{L}: number of confidence categories available for confidence ratings)
}

How the confidence variable \eqn{y} is computed varies across the different models.
The following models have been implemented so far:
\subsection{\strong{Signal detection rating model (SDT)}}{

According to SDT, the same sample of sensory
evidence is used to generate response and confidence, i.e.,
\eqn{y=x} and the confidence criteria span from the left and
right side of the decision criterion \eqn{c} (Green & Swets, 1966).
}

\subsection{\strong{Gaussian noise model (GN)}}{

According to the model, \eqn{y} is subject to
additive noise and assumed to be normally distributed around the decision
evidence value \eqn{x} with a standard deviation \eqn{\sigma} (Maniscalco & Lau, 2016).
The parameter  \eqn{\sigma} is a free parameter.
}

\subsection{\strong{Weighted evidence and visibility model (WEV)}}{

WEV assumes that the observer combines evidence about decision-relevant features
of the stimulus with the strength of evidence about choice-irrelevant features
to generate confidence (Rausch et al., 2018). Here, we use the version of the WEV model
used by Rausch et al. (2023), which assumes that \eqn{y} is normally
distributed with a mean of \eqn{(1-w)\times x+w \times d_k\times R} and standard deviation \eqn{\sigma}.
The parameter \eqn{\sigma} quantifies the amount of unsystematic variability
contributing to confidence judgments but not to the discrimination judgments.
The parameter \eqn{w} represents the weight that is put on the choice-irrelevant
features in the confidence judgment. \eqn{w} and \eqn{\sigma} are fitted in
addition to the set of shared parameters.
}

\subsection{\strong{Post-decisional accumulation model (PDA)}}{

PDA represents the idea of on-going information accumulation after the
discrimination choice (Rausch et al., 2018). The parameter \eqn{b} indicates the amount of additional
accumulation. The confidence variable is normally distributed with mean
\eqn{x+S\times d_k\times b} and variance \eqn{b}.
For this model the parameter \eqn{b} is fitted in
addition to the set of shared parameters.
}

\subsection{\strong{Independent Gaussian model (IG)}}{

According to IG, \eqn{y} is sampled independently
from \eqn{x} (Rausch & Zehetleitner, 2017). \eqn{y} is normally distributed with a mean of \eqn{a\times d_k} and variance
of 1 (again as it would scale with \eqn{m}). The free parameter \eqn{m}
represents the amount of information available for confidence judgment
relative to amount of evidence available for the discrimination decision and can
be smaller as well as greater than 1.
}

\subsection{\strong{Independent truncated Gaussian model: HMetad-Version (ITGc)}}{

According to the version of ITG consistent
with the HMetad-method (Fleming, 2017; see Rausch et al., 2023), \eqn{y} is sampled independently
from \eqn{x} from a truncated Gaussian distribution with a location parameter
of \eqn{S\times d_k \times m/2} and a scale parameter of 1. The Gaussian distribution of \eqn{y}
is truncated in a way that it is impossible to sample evidence that contradicts
the original decision: If \eqn{R = -1}, the distribution is truncated to the
right of \eqn{c}. If \eqn{R = 1}, the distribution is truncated to the left
of \eqn{c}. The additional parameter \eqn{m} represents metacognitive efficiency,
i.e., the amount of information available for confidence judgments relative to
amount of evidence available for discrimination decisions and  can be smaller
as well as greater than 1.
}

\subsection{\strong{Independent truncated Gaussian model: Meta-d'-Version (ITGcm)}}{

According to the version of the ITG consistent
with the original meta-d' method (Maniscalco & Lau, 2012, 2014; see Rausch et al., 2023),
\eqn{y} is sampled independently from \eqn{x} from a truncated Gaussian distribution with a location parameter
of \eqn{S\times d_k \times m/2} and a scale parameter
of 1. If \eqn{R = -1}, the distribution is truncated to the right of \eqn{m\times c}.
If \eqn{R = 1}, the distribution is truncated to the left of  \eqn{m\times c}.
The additional parameter \eqn{m} represents metacognitive efficiency, i.e.,
the amount of information available for confidence judgments relative to
amount of evidence available for the discrimination decision and  can be smaller
as well as greater than 1.
}

\subsection{\strong{Logistic noise model (logN)}}{

According to logN, the same sample
of sensory evidence is used to generate response and confidence, i.e.,
\eqn{y=x} just as in SDT (Shekhar & Rahnev, 2021). However, according to logN, the confidence criteria
are not assumed to be constant, but instead they are affected by noise drawn from
a lognormal distribution. In each trial, \eqn{\theta_{-1,i}} is given
by \eqn{c -  \epsilon_i}. Likewise,  \eqn{\theta_{1,i}} is given by
\eqn{c + \epsilon_i}. \eqn{\epsilon_i} is drawn from a lognormal distribution with
the location parameter
\eqn{\mu_{R,i}=log(|\overline{\theta}_{R,i}- c|) - 0.5 \times \sigma^{2}} and
scale parameter \eqn{\sigma}. \eqn{\sigma} is a free parameter designed to
quantify metacognitive ability. It is assumed that the criterion noise is perfectly
correlated across confidence criteria, ensuring that the confidence criteria
are always perfectly ordered. Because \eqn{\theta_{-1,1}}, ..., \eqn{\theta_{-1,L-1}},
\eqn{\theta_{1,1}}, ..., \eqn{\theta_{1,L-1}} change from trial to trial, they are not estimated
as free parameters. Instead, we estimate the means of the confidence criteria, i.e., \eqn{\overline{\theta}_{-1,1}, ...,
\overline{\theta}_{-1,L-1}, \overline{\theta}_{1,1}, ...  \overline{\theta}_{1,L-1}},
as free parameters.
}

\subsection{\strong{Logistic weighted evidence and visibility model (logWEV)}}{

logWEV is a combination of logN and WEV proposed by Shekhar and Rahnev (2023).
Conceptually, logWEV assumes that the observer combines evidence about decision-relevant features
of the stimulus with the strength of evidence about choice-irrelevant features (Rausch et al., 2018).
The model also assumes that noise affecting the confidence decision variable is lognormal
in accordance with Shekhar and Rahnev (2021).
According to logWEV, the confidence decision variable \eqn{y} is equal to
\eqn{y^*\times R}. \eqn{y^*} is sampled from a lognormal distribution with a location parameter
of \eqn{(1-w)\times x\times R + w \times d_k} and a scale parameter of \eqn{\sigma}.
The parameter \eqn{\sigma} quantifies the amount of unsystematic variability
contributing to confidence judgments but not to the discrimination judgments.
The parameter \eqn{w} represents the weight that is put on the choice-irrelevant
features in the confidence judgment. \eqn{w} and \eqn{\sigma} are fitted in
addition to the set of shared parameters.
}

}
}
\examples{
# 1. Select one subject from the masked orientation discrimination experiment
data <- subset(MaskOri, participant == 1)
head(data)

# 2. Use fitting function
\donttest{
  # Fitting takes some time (about 10 minutes on an 2.8GHz processor) to run:
  FitFirstSbjWEV <- fitConf(data, model="WEV")
}
}
\references{
Akaike, H. (1974). A New Look at the Statistical Model Identification. IEEE Transactions on Automatic Control, AC-19(6), 716–723.doi: 10.1007/978-1-4612-1694-0_16\cr

Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel inference: A practical information-theoretic approach. Springer.\cr

Fleming, S. M. (2017). HMeta-d: Hierarchical Bayesian estimation of metacognitive efficiency from confidence ratings. Neuroscience of Consciousness, 1, 1–14. doi: 10.1093/nc/nix007\cr

Green, D. M., & Swets, J. A. (1966). Signal detection theory and psychophysics. Wiley.\cr

Maniscalco, B., & Lau, H. (2012). A signal detection theoretic method for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.\cr

Maniscalco, B., & Lau, H. C. (2014). Signal Detection Theory Analysis of Type 1 and Type 2 Data: Meta-d’, Response- Specific Meta-d’, and the Unequal Variance SDT Model. In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp. 25–66). Springer. doi: 10.1007/978-3-642-45190-4_3\cr

Maniscalco, B., & Lau, H. (2016). The signal processing architecture underlying subjective reports of sensory awareness. Neuroscience of Consciousness, 1, 1–17. doi: 10.1093/nc/niw002\cr

Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in masked orientation judgments is informed by both evidence and visibility. Attention, Perception, and Psychophysics, 80(1), 134–154. doi: 10.3758/s13414-017-1431-5\cr

Rausch, M., Hellmann, S., & Zehetleitner, M. (2023). Measures of metacognitive efficiency across cognitive models of decision confidence. Psychological Methods. doi: 10.31234/osf.io/kdz34\cr

Rausch, M., & Zehetleitner, M. (2017). Should metacognition be measured by logistic regression? Consciousness and Cognition, 49, 291–312. doi: 10.1016/j.concog.2017.02.007\cr

Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461–464. doi: 10.1214/aos/1176344136\cr

Shekhar, M., & Rahnev, D. (2021). The Nature of Metacognitive Inefficiency in Perceptual Decision Making. Psychological Review, 128(1), 45–70. doi: 10.1037/rev0000249\cr

Shekhar, M., & Rahnev, D. (2023). How Do Humans Give Confidence? A Comprehensive Comparison of Process Models of Perceptual Metacognition. Journal of Experimental Psychology: General. doi:10.1037/xge0001524\cr

Peters, M. A. K., Thesen, T., Ko, Y. D., Maniscalco, B., Carlson, C., Davidson, M., Doyle, W., Kuzniecky, R., Devinsky, O., Halgren, E., & Lau, H. (2017). Perceptual confidence neglects decision-incongruent evidence in the brain. Nature Human Behaviour, 1(0139), 1–21. doi:10.1038/s41562-017-0139
}
\author{
Sebastian Hellmann, \email{sebastian.hellmann@tum.de}\cr
Manuel Rausch, \email{manuel.rausch@ku.de}
}
