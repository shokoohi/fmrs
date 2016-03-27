## ------------------------------------------------------------------------
library(fmrs)
set.seed(1980)
nComp = 2
nCov = 10
n = 500
REP = 500
deviance = c(1, 1)
pi = c(0.4, 0.6)
rho = 0.5
coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
umax = 40

## ------------------------------------------------------------------------
dat <- fmrs.gen.data(n = n, nComp = nComp, nCov = nCov,
                     coeff = c(coeff1, coeff2), deviance = deviance,
                     pi = pi, rho = rho, umax = umax, disFamily = "lnorm")

## ------------------------------------------------------------------------
res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
                    nComp = nComp, disFamily = "lnorm",
                    initCoeff = rnorm(nComp*nCov+nComp),
                    initDeviance = rep(1, nComp),
                    initPi = rep(1/nComp, nComp))

res.mle$coefficients
res.mle$deviance
res.mle$pi

## ------------------------------------------------------------------------
res.lam <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
                       nComp = nComp, disFamily = "lnorm",
                       initCoeff=c(res.mle$coefficients),
                       initDeviance = res.mle$deviance,
                       initPi = res.mle$pi, penFamily = "adplasso")
res.lam

## ------------------------------------------------------------------------
res.var <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
                       nComp = nComp, disFamily = "lnorm",
                       initCoeff = c(res.mle$coefficients),
                       initDeviance = res.mle$deviance,
                       initPi = res.mle$pi, penFamily = "adplasso",
                       lambPen = res.lam$lamPen)
res.var$coefficients
res.var$deviance
res.var$pi

## ------------------------------------------------------------------------
beta.est <- coefficients(res.var)[-1,]
round(beta.est,5)

