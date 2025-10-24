<h1 align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Optimal Nuisance Function Tuning</h1>
<h3 align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Optimal Nuisance Function Tuning for Estimating a Doubly Robust Functional under Proportional Asymptotics</h3>
<p align="center" style="margin-bottom:0px; border-bottom:0px; padding-bottom:0px">Sean McGrath, Debarghya Mukherjee, Rajarshi Mukherjee, Zixiao Jolene Wang</p>

<p align="center">
    <a style="text-decoration:none !important;" href="https://arxiv.org/abs/2509.25536" alt="arXiv"><img src="https://img.shields.io/badge/paper-arXiv-red" /></a>
</p>

We investigate asymptotically optimal tuning parameter selection in ridge regression for estimating nuisance functions within the **Expected Conditional Covariance (ECC)** framework.

Key contributions of our study include:

1. **Debiased ECC Estimators**:
   We derive debiased versions of three existing ECC estimators (Integral-based estimator, Newey Robins estimator, and doubly robust estimator) under linear models for (Y) and (A) on high-dimensional covariates (X \in \mathbb{R}^p), assuming the proportional asymptotic regime (p/n \to c \in (0, \infty)).

2. **√n-consistency via Ridge-based Nuisance Estimation**:
   By employing ridge regression to estimate nuisance functions and applying bias correction, our estimators achieve **√n-consistency** across multiple sample-splitting strategies.

3. **Asymptotic Variance Characterization**:
   We derive the asymptotic variances of the debiased estimators, highlighting the **interplay between sample splitting, estimator design, and tuning parameter choices**. Importantly, we show that prediction-optimal tuning parameters for nuisance estimation do **not necessarily minimize** the asymptotic variance of the ECC estimator.

We implemented two simulation settings:

#### Simulation 1: Asymptotic bias & variance of the integral-based estimator

1) We illustrate the asymptotic bias of the integral-based estimator, Newey Robins estimator, and doubly robust estimator and we verify that the debiased versions of these estimators are asymptotically unbiased.

2) We illustrate the asymptotic variances of the debiased estimators across a range of values for the nuisance parameters when $ c = 0.5$. We then compare the optimal tuning parameters to those that minimize the prediction error (i.e., the prediction-optimal tuning parameters), as in the simulations presented in the main text. For the sake of completeness, we also present the results for $c = 2$ where we do not trim the y-axis (unlike the main text).

3) We verify the derivations in Appendix F regarding the asymptotic mean squared error of the nuisance functions and the prediction-optimal tuning parameters.

The code is located in the folder  `exp1`, we generate Figure 1,2,3,4,5,6,7,8,9,10,11 based on this. 

#### Simulation 2. Estimating Logistic Model Coefficients and Quadratic Forms

We illustrate that the parametric bootstrap approach often performs well for estimating the asymptotic variance of the debiased estimators. 

The code is located in the folder  `exp2`, we generate Figure 12,13,14,15 based on this.

Any questions about the code, feel free to contact me at zixiaowang@fas.harvard.edu.


