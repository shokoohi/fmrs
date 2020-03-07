---
title: "Using fmrs package"
author: "Farhad Shokoohi"
date: "2020-02-23"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Sending Messages With Gmailr}
  %\usepackage[utf8]{inputenc}
---


# fmrs package in action

## Data generation
The function `fmrs.gendata` generates a data set from an FMRs. 
It has the form

````
fmrs.gendata(nObs, nComp, nCov, coeff, dispersion, mixProp, rho, umax, ...)
````

where `n` is the sample size, `nComp` is the order of FMRs, `nCov` is 
the number of regression covariates, `coeff`, `dispersion` and `mixProp` are 
the parameters of regression models, i.e. regression coefficients, 
dispersion (standard deviation) of the errors (sub-distributions) and 
mixing proportions, 
respectively, and `rho` is the used in the variance-covariance matrix for 
simulating the design matrix `x`, and `umax` is the upper bound for 
Uniform distribution for generating censoring times. 

Depending on the choice of `disFamily`, the function `fmrs.gendata` generates 
a simulated data from FMRs . For instance, if we choose 
`disFamily = "norm"`, the function ignores the censoring parameter `umax` 
and generates a data set from an FMR model with Normal sub-distributions. 
However, if we choose `disFamily = "lnorm"` or `disFamily = "weibull"`, the 
function generates data under a finite mixture of AFT regression model with 
Log-Normal or Weibull sub-distributions. 

The `fmrs.gendata` function returns a list which includes a vector of 
responses `$y`, a matrix of covariates `$x` and a vector of censoring 
indicators `$delta` as well as the name of sub-distribution of the 
mixture model.


## MLE of FMRs

The function `fmrs.mle` in fmrs package provides maximum likelihood 
estimation for the parameters of an FMRs. The function has the 
following form,

````
fmrs.mle(y, x, delta, nComp, ...)
````

where `y`, `x` and `delta` are observations, covariates and censoring 
indicators respectively, and `nComp` is the order of FMRs, `initCoeff`, 
`initDispersion` and `initmixProp` are initial values for EM and NR 
algorithms, and the rest of arguments of the function are controlling 
parameres. The function returns a fitted FMRs with estimates of 
regression parameters, standard deviations and mixing proportions of the 
mixture model.
It also returns the log-likelihood and BIC under the estimated model, the 
maximum number of iterations used in EM algorithm and the type of the 
fitted model.

Note that one can do ridge regression by setting a value for tuning 
parameter of the ridge penalty  other than zero in the argument `lambRidge`. 


## Variable selection in FMRs 

To do the variable selection we provided the function `fmrs.varsel` 
with the form

````
fmrs.varsel(y, x, delta, nComp, ...)
````

where `penFamily` is the penalty including `"adplasso"`, `"lasso"`, `"mcp"`, 
`"scad"`, `"sica"` and `"hard"`, and `lambPen` is the set of tuning 
parameters for components of penalty. We can run the function `fmrslme` 
first and use the parameter estimates as initial values for the 
function `fmrs.varsel`. 

### Choice of tuning parameter

There are two approaches for specifying tuning parameters: Common and 
Component-wise tuning parameters. If we consider choosing common tuning 
parameter, we can use the BIC criteria to search through the a set of 
candidate values in the interval $(0,\lambda_M)$, where $\lambda_M$ is a 
prespecified number. The BIC is provided by the function `fmrs.varsel` for 
each candidate point and we choose the optimal $\hat\lambda$, the one that 
minimizes BIC. This approach will be good for the situations with enough 
samples sizes. It also reduces the computational burden. 

On the other hand, if we consider choosing component-wise tuning parameters 
we use the following function to search for optimal 
$(\lambda_1, \ldots, \lambda_K)$ from the set of candidate values in 
$(0, \lambda_M)$. The function is

````
fmrs.tunsel(y, x, delta, nComp, ...)
````

It is a good practice run the function `fmrs.mle` first and use the parameter 
estimates as initial values in the function `fmrs.tunsel`.
The function `fmrs.mle` then returns a set optimal tuning parameters of 
size `nComp` to be used in variable selection function. Note that this 
approach still is under theoretical study and is not proved to give optimal 
values in general. 


## Example: finite mixture of AFT regression model (Log-Normal)

We use a simulated data set to illustrate using `fmrs` package. We generate 
the covariates (design matrix) from a multivariate normal distribution of 
dimension `nCov=10` and sample size 200 with mean vector $\bf 0$ and 
variance-covariance matrix $\Sigma=(0.5^{|i-j|})$. We then generate 
time-to-event data from a finite mixture of two components AFT regression 
models with Log-Normal sub-distributions. The mixing proportions are set 
to $\pi=(0.3, 0.7)$. We choose $\boldsymbol\beta_0=(2,-1)$ as the 
intercepts, $\boldsymbol\beta_1=(-1, -2, 1, 2, 0 , 0, 0, 0, 0, 0)$ and 
$\boldsymbol\beta_2=(1, 2, 0, 0, 0 , 0, -1, 2, -2, 3)$ as the regression 
coefficients in first and second component, respectively. 

We start by loading necessary libraries and assigning the parameters of 
model as follows.


```r
library(fmrs)
set.seed(1980)
nComp = 2
nCov = 10
nObs = 500
dispersion = c(1, 1)
mixProp = c(0.4, 0.6)
rho = 0.5
coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
umax = 40
```

Using the function `fmrs.gendata`, we generate a data set from a finite 
mixture of accelerated failure time regression models with above settings as 
follow. 


```r
dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
                     coeff = c(coeff1, coeff2), dispersion = dispersion,
                     mixProp = mixProp, rho = rho, umax = umax,
                     disFamily = "lnorm")
```

Now we assume that the generated data are actually real data. We find MLE of 
the parameters of the assumed model using following code. Note that almost 
all mixture of regression models depends on initial values. Here we 
generate the initial values form a Normal distribution with 


```r
res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
                   nComp = nComp, disFamily = "lnorm",
                   initCoeff = rnorm(nComp*nCov+nComp),
                   initDispersion = rep(1, nComp),
                   initmixProp = rep(1/nComp, nComp))
coefficients(res.mle)
```

```
##                Comp.1     Comp.2
## Intercept -0.26344394  0.9976953
## X.1        0.32900425 -1.2803480
## X.2        0.06690871  0.1928663
## X.3        0.69669344  0.7978279
## X.4        0.51597721  0.4171002
## X.5       -0.09399042  0.1369026
## X.6        0.26679989 -0.1871274
## X.7       -0.17686523 -0.7404996
## X.8       -0.63562498  1.2610768
## X.9        1.23776252  1.1937360
## X.10      -2.20704095 -0.1349468
```

```r
dispersion(res.mle)
```

```
##        Comp.1   Comp.2
## [1,] 2.572907 3.213501
```

```r
mixProp(res.mle)
```

```
##         Comp.1    Comp.2
## [1,] 0.5693449 0.4306551
```

As we see the ML estimates of regression coefficients are not zero. 
Therefore MLE cannot provide a sparse solution. In order to provide a sparse 
solution, we use the variable selection methods developed by 
Shokoohi et. ale. (2016). First we need to find a good set of tuning 
parameters.
It can be done by using component-wise tuning parameter selection function 
implemented in `fmrs` as follows. In some settings, however, it is better to 
investigate if common tuning parameter performs better. 


```r
res.lam <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
                      nComp = nComp, disFamily = "lnorm",
                      initCoeff = c(coefficients(res.mle)),
                      initDispersion = dispersion(res.mle),
                      initmixProp = mixProp(res.mle),
                      penFamily = "adplasso")
show(res.lam)
```

```
## An object of class 'fmrstunpar'
##   Finite Mixture of Accelerated Failure Time Regression Models
##     Component-Wise Distributions: Log-Normal 
##   No. Components: 2
##   Penalty: adplasso
```

We have used MLE estimates as initial values for estimating the tuning 
parameters. Now we used the same set of values to do variable selection with 
adaptive lasso penalty as follows. 


```r
res.var <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
                      nComp = ncomp(res.mle), disFamily = "lnorm",
                      initCoeff=c(coefficients(res.mle)),
                      initDispersion = dispersion(res.mle),
                      initmixProp = mixProp(res.mle),
                      penFamily = "adplasso",
                      lambPen = slot(res.lam, "lambPen"))
coefficients(res.var)
```

```
##                  Comp.1        Comp.2
## Intercept -1.117163e+00  1.961682e+00
## X.1       -7.524812e-01  1.195135e+00
## X.2        3.975668e-01 -4.163286e-13
## X.3        2.218990e+00 -5.756995e-01
## X.4        2.618467e-02  1.650138e-01
## X.5        4.984195e-12  9.354113e-12
## X.6       -5.546534e-12  3.799180e-12
## X.7       -5.055157e-11  2.415133e-11
## X.8       -9.047182e-01  5.849482e-01
## X.9        1.951799e+00  6.691085e-03
## X.10      -2.083846e+00  4.842384e-14
```

```r
dispersion(res.var)
```

```
##         Comp.1   Comp.2
## [1,] 0.9021074 2.568686
```

```r
mixProp(res.var)
```

```
##         Comp.1    Comp.2
## [1,] 0.5536024 0.4463976
```

The final variables that are selected using this procedure are those with 
non-zero coefficients. Note that a switching between components of mixture 
has happened here.   


```r
round(coefficients(res.var)[-1,],3)
```

```
##      Comp.1 Comp.2
## X.1  -0.752  1.195
## X.2   0.398  0.000
## X.3   2.219 -0.576
## X.4   0.026  0.165
## X.5   0.000  0.000
## X.6   0.000  0.000
## X.7   0.000  0.000
## X.8  -0.905  0.585
## X.9   1.952  0.007
## X.10 -2.084  0.000
```

Therefore the variable selection and parameter estimation is done 
simultaneously using the fmrs package. 




