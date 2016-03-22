#' @title Variable Selection in Finite Mixture of Accelerated Failure Time Regression Models
#'  and Finite Mixture of Regression Models
#'
#' @description It provides variable selection and parameter estimation for Finite Mixture of Regression Models
#' and Finite Mixture of Accelerated Failure Time Regression Models with Log-Normal or Weibull sub-distributions.
#' The penalties lasso, scad, mcp, sica, adaptive lasso and hard are implemented. It also provide Ridge Regression and Elastic Net.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.varsel
#' @param y Time-to-event observations
#' @param x Design matrix (covariates)
#' @param delta Censoring indicator vector
#' @param nComp Orde (Number of components) of mixture model
#' @param disFamily Specify sub-distributions family. The options are \code{"norm"} for FMR models,
#' \code{"lnorm"} for mixture of AFT regression models with Log-Normal sub-distributions,
#' \code{"weibull"} for mixture of AFT regression models with Weibull sub-distributions,
#' @param initCoeff Vector of initial values for coefficients including intercepts.
#' @param initSigma Vector of initial values for standard deviations
#' @param initPi Vector of initial values for proportion of components
#' @param penFamily The penalty that used in variable selection method.
#' The available options are  \code{"lasso"}, \code{"adplasso"}, \code{"mcp"}, \code{"scad"}, \code{"sica"} and \code{"hard"}.
#' @param lambPen A vector of lambda for penaly
#' @param lambRidge Lambda for ridge penalty
#' @param nIterEM Number of iteration for EM algorithm
#' @param nIterNR Number of iteration for Newton-Raphson algorithm
#' @param conveps A positive number for avoiding NaN in computing divisions.
#' @param convepsEM Treshold for convergence of EM algorithm
#' @param convepsNR Treshold for convergence of Newton-Raphson algorithm
#' @param porNR Used in pow(0.5, porNR) for tuning the increment in Newton-Raphson algorithm.
#' @param gamMixPor Proportion of mixing parameters in the penalty. The value must belong to the interval [0,1]. If \code{gamMixPor = 0}, the penalty structure is no longer mixture.
#' @keywords FMR, AFT, Censored Data, EM Algorithm, Ridge Regression
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S. (2016) Variable Selection in Mixture of Survival Models
#' @return Parameter estimates of mixture of Normal or mixture of Log-Normal AFT models
#' @return A fitted FMRs model object of class \code{\link{fmrs.fit}}
#' @examples \dontrun{Variable Selection in Ovarian Cancer analysis
#'
#' }
#' @export
fmrs.varsel <- function(y,
                        x,
                        delta,
                        nComp,
                        disFamily = "lnorm",
                        initCoeff,
                        initSigma,
                        initPi,
                        penFamily = "lasso",
                        lambPen,
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

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("Cov",1:nCov,sep=".")))
  } else{
    xnames <- c("Intercept",colnames(x))
  }
  comnames <- c(paste("Comp",1:nComp,sep="."))

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

    res=.C("FMR_Norm_Surv_EM_VarSel", PACKAGE="fmrs",
           y = as.double(y),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Pen = as.double(lambPen),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           NumIterationEM = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
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
           GIC = as.double(0)
    )

  }else if(disFamily == "lnorm"){
    meth = "FMAFTR"
    logy = log(y)

    res=.C("FMR_Norm_Surv_EM_VarSel", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Pen = as.double(lambPen),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           NumIterationEM = as.integer(nIterEM),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
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
           GIC = as.double(0)
    )
  }else if(disFamily == "weibull"){
    meth = "FMAFTR"
    logy = log(y)

    res=.C("FMR_Weibl_Surv_EM_VarSel", PACKAGE="fmrs",
           y = as.double(logy),
           x = as.double(as.vector(unlist(x))),
           delta = as.double(delta),
           myPenaltyFamily = as.integer(myPenaltyFamily),
           Lambda.Pen = as.double(lambPen),
           Lambda.Ridge = as.double(lambRidge),
           Num.Comp = as.integer(nComp),
           Num.Cov = as.integer(nCov),
           Sample.Size = as.integer(n),
           NumIterationEM = as.integer(nIterEM),
           NumIterationNR = as.integer(nIterNR),
           PortionNF = as.integer(porNR),
           Max.iterEM.used = as.integer(0),
           Initial.Intercept = as.double(c(coef0[,1])),
           Initial.Coefficient = as.double(c(t(coef0[,-1]))),
           Initial.Sigma = as.double(initSigma),
           Initial.Pi = as.double(initPi),
           conv.eps = as.double(convepsNR),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
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
           GIC = as.double(0)
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
              penFamily = penFamily,
              lamPen = lambPen,
              lambRidge = lambRidge
  )
  class(fit) <- "fmrs.fit"
  return(fit)
}
