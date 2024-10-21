### Information-theoretic measures of metacognition

Dayan (2023) proposed several measures of metacognition based on
quantities of information theory.

- Meta-I is a measure of metacognitive sensitivity defined as the mutual
  information between the confidence and accuracy and is calculated as
  the transmitted information minus the minimal information given the
  accuracy:

$$meta-I = I(Y; \hat{Y}, C) - I(Y; \hat{Y})$$ This is equivalent to
Dayan’s formulation where meta-I is the information that confidences
transmit about the correctness of a response:

$$meta-I = I(Y = \hat{Y}; C)$$ - Meta-$I_{1}^{r}$ is meta-I normalized
by the value of meta-I expected assuming a signal detection model (Green
& Swets, 1966) with Gaussian noise, based on calculating the sensitivity
index d’:

$$meta-I_{1}^{r} = meta-I / meta-I(d')$$ - Meta-$I_{2}^{r}$ is meta-I
normalized by its theoretical upper bound, which is the information
entropy of accuracy, $H(Y = \hat{Y})$:

$$meta-I_{2}^{r} = meta-I / H(Y = \hat{Y})$$

Notably, Dayan (2023) pointed out that a liberal or conservative use of
the confidence levels will affected the mutual information and thus all
information-theoretic measures of metacognition.

In addition to Dayan’s measures, Meyen et al. (submitted) suggested an
additional measure that normalizes the Meta-I by the range of possible
values it can take. This required deriving lower and upper bound of the
transmitted information given a participant’s accuracy.

$$RMI = \frac{meta-I}{\max_{\text{accuracy}}\{meta-I\}}$$

As these measures are prone to estimation bias, the package offers a
simple bias reduction mechanism in which the observed frequencies of
stimulus-response combinations are taken as the underyling probability
distribution. From this, Monte-Carlo simulations are conducted to
estimate and subtract the bias in these measures. Note that there
provably is no way to completely remove this bias.
