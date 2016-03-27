#' @title  Component-Wise Tuning Parameter Selection in Finite Mixture of
#'     Accelerated Failure Time Regression Models
#'     and Finite Mixture of Regression Models
#'
#' @description  It provides component-wise tuning parameters for Finite
#'     Mixture of Accelerated Failure Time Regression Models
#'     and Finite Mixture of Regression Models.
#'     The penalties that are implemented in this package are \code{lasso},
#'     \code{adplasso}, \code{scad}, \code{mcp}, \code{sica} and \code{hard}.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrstunsel
#' @param y Responses (observations)
#' @param x Design matrix (covariates)
#' @param delta Censoring indicator vector
#' @param nComp Order (Number of components) of mixture model
#' @param disFamily Specify sub-distributions family. The options
#'     are \code{"norm"} for FMR models,
#'     \code{"lnorm"} for mixture of AFT regression models with Log-Normal
#'     sub-distributions, \code{"weibull"} for mixture of AFT regression models
#'     with Weibull sub-distributions,
#' @param initCoeff Vector of initial values for regression coefficients
#'     including intercepts
#' @param initDeviance Vector of initial values for standard deviations
#' @param initPi Vector of initial values for proportion of components
#' @param penFamily Penalty name that is used in variable selection method.
#'     The available options are  \code{"lasso"}, \code{"adplasso"},
#'     \code{"mcp"}, \code{"scad"}, \code{"sica"} and \code{"hard"}.
#' @param lambRidge A positive value for tuniing parameter in Ridge regression
#'     or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in
#'     EM algorithm
#' @param convepsNR A positive value for treshold of convergence in
#'     NR algorithm
#' @param porNR Used in pow(0.5, porNR) for tuning the increment in
#'     NR algorithm
#' @param gamMixPor Proportion of mixing parameters in the penalty. The value
#'     must be in the interval [0,1]. If \code{gamMixPor = 0}, the penalty
#'     structure is no longer mixture.
#' @keywords FMR, AFT, Censored Data, EM Algorithm, Ridge Regression
#' @concept fmr, aft, lasso, adplasso, mcp, scad, sica, ridge
#' @details The maximizer of penalized Log-Likelihood depends on selecting a
#'     set of good tuning parameters which is a rather thorny issue. We choose
#'     a value in an equally spaced set of values in \eqn{(0, \lambda_{max})}
#'     for a pre-specified \eqn{\lambda_{max}} that maximize the component-wise
#'     BIC, \deqn{\hat\lambda_{k} ={argmax}_{\lambda_{k}}BIC_k(\lambda_{k}) =
#'     {argmax}_{\lambda_{k}}\left\{\ell^{c}_{k, n}
#'     (\hat{\boldsymbol\Psi}_{\lambda_{k}, k}) -
#'     |d_{\lambda_{k},k}| \log (n)\right\},}
#'     where \eqn{d_{\lambda_{k},k}=\{j:\hat{\beta}_{\lambda_{k},kj}\neq 0,
#'     j=1,\ldots,d\}} is the active set  excluding the intercept
#'     and \eqn{|d_{\lambda_{k},k}|}
#'     is its size. This approach is much faster than using an \code{nComp}
#'     by \code{nComp} grid to select the set \eqn{\boldsymbol\lambda} to
#'     maximize the penallized Log-Likelihood.
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S.
#'     (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{tunepar-class}} object includes component-wise
#'     tuning parameter estimates to be used in variable selection procedure.
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
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      pi = pi, rho = rho, umax = umax, disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                     nComp = nComp, disFamily = "lnorm",
#'                     initCoeff = rnorm(nComp*nCov+nComp),
#'                     initDeviance = rep(1, nComp),
#'                     initPi = rep(1/nComp, nComp))
#'
#' res.lam <- fmrstunsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                        nComp = nComp, disFamily = "lnorm",
#'                        initCoeff=c(res.mle$coefficients),
#'                        initDeviance = res.mle$deviance,
#'                        initPi = res.mle$pi, penFamily = "adplasso")
#'
#' res.lam
#' @export
fmrstunsel <- function(y,
                        x,
                        delta,
                        nComp,
                        disFamily = "lnorm",
                        initCoeff,
                        initDeviance,
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
  if(is.null(initCoeff) | is.null(initDeviance) | is.null(initPi))
    stop("Initial values are not specified.")
  if(length(initCoeff) != nComp*nCov+nComp | length(initPi)!=nComp |
     length(initDeviance)!=nComp)
    stop("The length of initial values are not correctly specified.")
  if(!is.matrix(x))
    stop("Provide a matix for covariates.")
  nCov = dim(x)[2]
  n = length(y)
  if(dim(x)[1]!=n)
    stop("The length of observations and rows of design
         matrix does not match.")

  coef0 <- matrix(initCoeff, nrow = nComp, ncol = nCov+1, byrow = TRUE)

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
           Initial.Deviance = as.double(initDeviance),
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
           Initial.Deviance = as.double(initDeviance),
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
           Initial.Deviance = as.double(initDeviance),
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


  lambdafit <- list(lamPen =
                      matrix(res$Opt.Lambda, nrow = 1,
                             dimnames = c(list(NULL,c(paste("Comp",1:nComp,sep
                                                            = "."))))),
                    disFamily = disFamily,
                    penFamily = penFamily,
                    lamRidge = lambRidge,
                    method = meth
  )
  class(lambdafit) <- "tunepar"
  return(lambdafit)
}
