#' @title  Maximum Likelihood Estimation in Finite Mixture of Accelerated Failure Time Regression Models
#' and Finite Mixture of Regression Models
#'
#' @description  It provides parameter estimation for Finite Mixture of Accelerated Failure Time Regression Models and Finite Mixture of Regression Models.
#' It also provide Ridge Regression.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.mle
#' @param y Responses (observations)
#' @param x Design matrix (covariates)
#' @param delta Censoring indicator vector
#' @param nComp Order (Number of components) of mixture model
#' @param disFamily Specify sub-distributions family. The options are \code{"norm"} for FMR models,
#' \code{"lnorm"} for mixture of AFT regression models with Log-Normal sub-distributions,
#' \code{"weibull"} for mixture of AFT regression models with Weibull sub-distributions,
#' @param initCoeff Vector of initial values for regression coefficients including intercepts
#' @param initDeviance Vector of initial values for standard deviations
#' @param initPi Vector of initial values for proportion of components
#' @param lambRidge A positive value for Lambda in Ridge regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in EM algorithm
#' @param convepsNR A positive value for treshold of convergence in Newton-Raphson algorithm
#' @param porNR Used in pow(0.5, porNR) for tuning the increment in Newton-Raphson algorithm.
#' @keywords FMR, AFT, Censored Data, EM Algorithm, Ridge Regression
#' @concept fmr, aft, lasso, adplasso, mcp, scad, sica, ridge
#' @details Finite mixture of AFT regression models are represented as follows.
#' Let \eqn{X} be the survival time with non-negative values, and \eqn{\boldsymbol{z} =(z_{1}, \ldots, z_{d})^{\top}}
#' be a \eqn{d}-dimensional vector of covariates that may have an effect on \eqn{X}.
#' If the survival time is subject to right censoring, then the observed response time is \eqn{T=\min \{Y, C\}},
#' where \eqn{Y=\log X}, \eqn{C} is  logarithm of  the censoring time and \eqn{\delta=I_{\{y<c\}}} is the censoring indicator.
#' We say that \eqn{V=(T,\delta,\boldsymbol z)} follows a finite mixture of AFT regression models of order \eqn{K}
#' if the conditional density of \eqn{(T,\delta)} given \eqn{\boldsymbol z} has the form
#' \deqn{f(t,\delta;\boldsymbol{z},\boldsymbol\Psi)=\sum\limits_{k=1}^{K}\pi_{k}[f_Y(t;\theta_{k}(\boldsymbol z),
#' \sigma_{k})]^{\delta}[S_Y(t;\theta_{k}(\boldsymbol z),\sigma_{k})]^{1-\delta}[f_{C}(t)]^{1-\delta}[S_{C}(t)]^{\delta}}
#' where \eqn{f_Y(.)} and \eqn{S_Y(.)} are respectively the density and survival functions of \eqn{Y},
#' \eqn{f_C(.)} and \eqn{S_C(.)} are respectively the density and survival functions of \eqn{C};
#' and \eqn{{\theta}_{k}(\boldsymbol{z})=h(\beta_{0k}+\boldsymbol{z}^{\top}\boldsymbol\beta_{k})}
#' for a known link function \eqn{h(.)},  \eqn{\boldsymbol\Psi=(\pi_{1},\ldots,\pi_{K},\beta_{01},\ldots,
#' \beta_{0K}, \boldsymbol\beta_{1},\ldots,\boldsymbol\beta_{K},\sigma_{1},\ldots,\sigma_{K})^{\top}}
#' with \eqn{\boldsymbol\beta_{k}=(\beta_{k1},\beta_{k2},\ldots,\beta_{kd})^{\top}} and \eqn{0<\pi_{k}<1} with \eqn{\sum_{k=1}^{K}\pi_{k}=1}.
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S. (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{fmrs.fit-class}} object which includes parameter estimates of an FMRs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' REP = 500
#' deviance = c(1, 1)
#' pi = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gen.data(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      pi = pi, rho = rho, umax = umax, disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                     nComp = nComp, disFamily = "lnorm",
#'                     initCoeff = rnorm(nComp*nCov+nComp),
#'                     initDeviance = rep(1, nComp),
#'                     initPi = rep(1/nComp, nComp))
#'
#' res.mle$coefficients
#' res.mle$deviance
#' res.mle$pi
#' @export
fmrs.mle <- function(y,
                     x,
                     delta,
                     nComp,
                     disFamily = "lnorm",
                     initCoeff,
                     initDeviance,
                     initPi,
                     lambRidge = 0,
                     nIterEM = 2000,
                     nIterNR = 2,
                     conveps = 1e-8,
                     convepsEM = 1e-8,
                     convepsNR = 1e-8,
                     porNR = 2
)
{
  if(is.null(y))
    stop("Response vector is not specified.")
  if(is.null(x))
    stop("The desing matrix is not specified.")
  if(is.null(delta) & (disFamily!="norm"))
    stop("The censoring vector is not specified.")
  if(is.null(nComp))
    stop("Number of components of mixture model is not specified.")
  if(is.null(initCoeff) | is.null(initDeviance) | is.null(initPi))
    stop("Initial values are not specified.")
  if(length(initCoeff) != nComp*nCov+nComp | length(initPi)!=nComp | length(initDeviance)!=nComp)
    stop("The length of initial values are not correctly specified.")
  if(!is.matrix(x))
    stop("Provide a matix for covariates.")
  nCov = dim(x)[2]
  n = length(y)
  if(dim(x)[1]!=n)
    stop("The length of observations and rows of design matrix does not match.")

  coef0 <- matrix(initCoeff, nrow = nComp, ncol = nCov+1, byrow = T)

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("Cov",1:nCov,sep=".")))
  } else{
    xnames <- c("Intercept",colnames(x))
  }
  comnames <- c(paste("Comp",1:nComp,sep="."))

  if(disFamily == "norm"){
    meth = "FMR"
    delta = rep(1, n)
    res=.C("FMR_Norm_Surv_EM_MLE", PACKAGE="fmrs",
           y = as.double(y),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Num.iteration = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(unlist(coef0[,1])),
           Initial.Coefficient = as.double(unlist(t(coef0[,-1]))),
           Initial.Deviance = as.double(initDeviance),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Deviance.Hat = as.double(rep(0,nComp)),
           Pi.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )
  }else if(disFamily == "lnorm"){
    meth = "FMAFTR"
    logy = log(y)
    res=.C("FMR_Norm_Surv_EM_MLE", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Num.iteration = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(unlist(coef0[,1])),
           Initial.Coefficient = as.double(unlist(t(coef0[,-1]))),
           Initial.Deviance = as.double(initDeviance),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Deviance.Hat = as.double(rep(0,nComp)),
           Pi.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )

  }else if(disFamily == "weibull"){
    meth = "FMAFTR"
    logy = log(y)
    res=.C("FMR_Weibl_Surv_EM_MLE", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Num.iterationEM = as.integer(nIterEM),
           Num.iterationNR = as.integer(nIterNR),
           PortionNF = as.integer(porNR),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Deviance = as.double(initDeviance),
           Initial.Pi = as.double(initPi),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Deviance.Hat = as.double(rep(0,nComp)),
           Pi.Hat = as.double(rep(0,nComp)),
           LogLikelihood = as.double(0),
           BIC = as.double(0),
           AIC = as.double(0),
           GCV = as.double(0),
           EBIC1 = as.double(0),
           EBIC5 = as.double(0),
           GIC = as.double(0),
           predict = as.double(rep(0,n*nComp)),
           residual = as.double(rep(0,n*nComp)),
           tau = as.double(rep(0,n*nComp))
    )
  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }

  fit <- list(coefficients = array(rbind(res$Intecept.Hat, matrix(res$Coefficient.Hat, nrow = nCov, byrow = F)),
                                   dim = c(nCov+1, nComp), dimnames = list(xnames,comnames)),
              deviance = array(res$Deviance.Hat, dim = c(1,nComp),dimnames = list(NULL,comnames)),
              pi = array(res$Pi.Hat, dim = c(1,nComp),dimnames = list(NULL,comnames)),
              logLik = res$LogLikelihood,
              BIC = res$BIC,
              nIterEMconv = res$Max.iterEM.used,
              method = meth,
              disFamily = disFamily,
              lambRidge = lambRidge,
              fitted = array(matrix(res$predict, nrow = n, byrow = F),
                             dim = c(n, nComp), dimnames = list(NULL,comnames)),
              residuals = array(matrix(res$residual, nrow = n, byrow = F),
                                dim = c(n, nComp), dimnames = list(NULL,comnames)),
              weights = array(matrix(res$tau, nrow = n, byrow = F),
                              dim = c(n, nComp), dimnames = list(NULL,comnames))
  )
  class(fit) <- "fmrs.fit"
  return(fit)

}
