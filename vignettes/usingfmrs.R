## ------------------------------------------------------------------------
library(fmrs)
set.seed(1980)
nComp = 2
nCov = 10
n = 500
REP = 500
deviance = c(1, 1)
mixProp = c(0.4, 0.6)
rho = 0.5
coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
umax = 40

## ------------------------------------------------------------------------
dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
                     coeff = c(coeff1, coeff2), deviance = deviance,
                     mixProp = mixProp, rho = rho, umax = umax, 
                     disFamily = "lnorm")

## ------------------------------------------------------------------------
res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
                    nComp = nComp, disFamily = "lnorm",
                    initCoeff = rnorm(nComp*nCov+nComp),
                    initDeviance = rep(1, nComp),
                    initmixProp = rep(1/nComp, nComp))

coefficients(res.mle)
deviance(res.mle)
mixProp(res.mle)

## ------------------------------------------------------------------------
res.lam <- fmrstunsel(y = dat$y, x = dat$x, delta = dat$delta,
                       nComp = nComp, disFamily = "lnorm",
                       initCoeff = coefficients(res.mle),
                       initDeviance = deviance(res.mle),
                       initmixProp = mixProp(res.mle), 
                       penFamily = "adplasso")
slot(res.lam, "lambPen")

## ------------------------------------------------------------------------
res.var <- fmrsvarsel(y = dat$y, x = dat$x, delta = dat$delta,
                       nComp = nComp, disFamily = "lnorm",
                       initCoeff = coefficients(res.mle),
                       initDeviance = deviance(res.mle),
                       initmixProp = mixProp(res.mle), 
                       penFamily = "adplasso",
                       lambPen = slot(res.lam, "lambPen"))
coefficients(res.var)
deviance(res.var)
mixProp(res.var)

## ------------------------------------------------------------------------
round(coefficients(res.var)[-1,],5)

