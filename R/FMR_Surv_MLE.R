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
#' @param initSigma Vector of initial values for standard deviations
#' @param initPi Vector of initial values for proportion of components
#' @param lambRidge A positive value for Lambda in Ridge regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in EM algorithm
#' @param convepsNR A positive value for treshold of convergence in Newton-Raphson algorithm
#' @param porNR Used in pow(0.5, porNR) for tuning the increment in Newton-Raphson algorithm.
#' @keywords FMR, AFT, Censored Data, EM Algorithm, Ridge Regression
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S. (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{fmrs.fit}} object which includes parameter estimates of an FMRs model
#' @examples \dontrun{ MLE of the generated data, see fmrs.gen.data
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                     nComp = nComp, disFamily = "lnorm",
#'                     initCoeff = rnorm(nComp*nCov+nComp),
#'                     initSigma = rep(1, nComp),
#'                     initPi = rep(1/nComp, nComp))
#'
#' res.mle$coefficients
#' res.mle$sigma
#' res.mle$pi
#' }
#' @export
fmrs.mle <- function(y,
                     x,
                     delta,
                     nComp,
                     disFamily = "lnorm",
                     initCoeff,
                     initSigma,
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
  if(is.null(initCoeff) | is.null(initSigma) | is.null(initPi))
    stop("Initial values are not specified.")
  if(length(initCoeff) != nComp*nCov+nComp | length(initPi)!=nComp | length(initSigma)!=nComp)
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
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Sigma.Hat = as.double(rep(0,nComp)),
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
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Sigma.Hat = as.double(rep(0,nComp)),
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
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Sigma.Hat = as.double(rep(0,nComp)),
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
              sigma = array(res$Sigma.Hat, dim = c(1,nComp),dimnames = list(NULL,comnames)),
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
