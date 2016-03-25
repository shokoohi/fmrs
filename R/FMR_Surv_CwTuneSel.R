#' @title  Component-Wise Tuning Parameter Selection in Finite Mixture of Accelerated Failure Time Regression Models
#' and Finite Mixture of Regression Models
#'
#' @description  It provides component-wise tuning parameters for Finite Mixture of Accelerated Failure Time Regression Models
#' and Finite Mixture of Regression Models.
#' The penalties that are implemented in this package are \code{lasso}, \code{adplasso}, \code{scad}, \code{mcp}, \code{sica} and \code{hard}.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.tunsel
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
#' @param penFamily Name of the penalty that is used in variable selection method.
#' The available options are  \code{"lasso"}, \code{"adplasso"}, \code{"mcp"}, \code{"scad"}, \code{"sica"} and \code{"hard"}.
#' @param lambRidge A positive value for Lambda in Ridge regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in EM algorithm
#' @param convepsNR A positive value for treshold of convergence in Newton-Raphson algorithm
#' @param porNR Used in pow(0.5, porNR) for tuning the increment in Newton-Raphson algorithm.
#' @param gamMixPor Proportion of mixing parameters in the penalty. The value must be in the interval [0,1]. If \code{gamMixPor = 0}, the penalty structure is no longer mixture.
#' @keywords FMR, AFT, Censored Data, EM Algorithm, Ridge Regression
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S. (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{fmrs.lambda}} object which includes component-wise tuning parameter estimates to be used in variable selection procedure.
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' REP = 500
#' sigma = c(1, 1)
#' pi = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gen.data(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), sigma = sigma,
#'                      pi = pi, rho = rho, umax = umax, disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                     nComp = nComp, disFamily = "lnorm",
#'                     initCoeff = rnorm(nComp*nCov+nComp),
#'                     initSigma = rep(1, nComp),
#'                     initPi = rep(1/nComp, nComp))
#'
#' res.lam <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                        nComp = nComp, disFamily = "lnorm",
#'                        initCoeff=c(res.mle$coefficients),
#'                        initSigma = res.mle$sigma,
#'                        initPi = res.mle$pi, penFamily = "adplasso")
#'
#' res.lam
#' @export
fmrs.tunsel <- function(y,
                        x,
                        delta,
                        nComp,
                        disFamily = "lnorm",
                        initCoeff,
                        initSigma,
                        initPi,
                        penFamily = "lasso",
                        lambRidge = 0,
                        nIterEM = 2000,
                        nIterNR = 2,
                        conveps = 1e-8,
                        convepsEM = 1e-8,
                        convepsNR = 1e-8,
                        porNR = 2,
                        gamMixPor = 1
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

  if(penFamily == "lasso") myPenaltyFamily = 1
  else if (penFamily == "scad") myPenaltyFamily = 2
  else if (penFamily == "mcp") myPenaltyFamily = 3
  else if (penFamily == "sica") myPenaltyFamily = 4
  else if (penFamily == "adplasso") myPenaltyFamily = 5
  else if (penFamily == "hard") myPenaltyFamily = 6
  else {print("Penalty is not specified.")  }

  if(disFamily == "norm"){
    meth = "FMR"
    delta = rep(1, n)

    res=.C("FMR_Norm_Surv_CwTuneParSel", PACKAGE="fmrs",
           y = as.double(y),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp))
    )
  }else if(disFamily == "lnorm"){
    meth = "FMAFTR"
    logy = log(y)
    res=.C("FMR_Norm_Surv_CwTuneParSel", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp))
    )
  }else if(disFamily == "weibull"){
    meth = "FMAFTR"
    logy = log(y)

    res=.C("FMR_Weibl_Surv_CwTuneParSel", PACKAGE = "fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           Num.NRiteration = as.double(nIterNR),
           Num.PortionNF = as.double(porNR),
           conv.eps = as.double(convepsNR),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp))
    )
  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }


  lambdafit <- list(lamPen = matrix(res$Opt.Lambda, nrow = 1, dimnames = c(list(NULL,c(paste("Comp",1:nComp,sep = "."))))),
                    disFamily = disFamily,
                    penFamily = penFamily,
                    lamRidge = lambRidge,
                    method = meth
  )
  class(lambdafit) <- "fmrs.lambda"
  return(lambdafit)
}
