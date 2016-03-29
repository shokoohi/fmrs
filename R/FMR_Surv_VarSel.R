#' @title Variable Selection in Finite Mixture of Accelerated Failure Time
#'     Regression Models and Finite Mixture of Regression Models
#'
#' @description It provides variable selection and parameter estimation for
#'     Finite Mixture of Accelerated Failure Time Regression (FMAFTR) Models
#'     and Finite Mixture of Regression (FMR) Models.
#'     The penalties that are implemented in this package are \code{lasso},
#' \code{adplasso}, \code{scad}, \code{mcp}, \code{sica} and \code{hard}.
#'     It also provide Ridge Regression and Elastic Net.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrsvarsel
#' @param y Responses (observations)
#' @param x Design matrix (covariates)
#' @param delta Censoring indicators
#' @param nComp Order (Number of components) of mixture model
#' @param disFamily Name of sub-distributions' family. The options
#'     are \code{"norm"} for FMR models, \code{"lnorm"} for mixture of AFT
#'     regression models with Log-Normal sub-distributions, \code{"weibull"}
#'     for mixture of AFT regression models with Weibull sub-distributions
#' @param initCoeff Vector of initial values for regression coefficients
#'     including intercepts
#' @param initDeviance Vector of initial values for standard deviations
#' @param initmixProp Vector of initial values for proportion of components
#' @param penFamily Penalty name that is used in variable selection method
#'     The available options are  \code{"lasso"}, \code{"adplasso"},
#'     \code{"mcp"}, \code{"scad"}, \code{"sica"} and \code{"hard"}.
#' @param lambPen A vector of positive numbers for tuning parameters
#' @param lambRidge A positive value for tuning parameter in Ridge
#'     Regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in
#'     EM algorithm
#' @param convepsNR A positive value for treshold of convergence in
#'     NR algorithm
#' @param porNR Used in pow(0.5, porNR) for tuning the increment in
#'     NR algorithm
#' @param gamMixPor Proportion of mixing parameters in the penalty. The
#'     value must be in the interval [0,1]. If \code{gamMixPor = 0}, the
#'     penalty structure is no longer mixture.
#' @keywords FMR, AFT, Censored Data, EM Algorithm, Ridge Regression
#' @concept fmr, aft, lasso, adplasso, mcp, scad, sica, ridge
#' @details The penalized likelihood of a finite mixture of AFT regression
#'     models is written as \deqn{\tilde\ell_{n}(\boldsymbol\Psi)
#'     =\ell_{n}(\boldsymbol\Psi) -
#'     \mathbf{p}_{\boldsymbol\lambda_{n}}(\boldsymbol\Psi)}
#'     where \deqn{\mathbf{p}_{\boldsymbol\lambda_{n}}(\boldsymbol\Psi) =
#'     \sum\limits_{k=1}^{K}\pi_{k}^\alpha\left\{
#'     \sum\limits_{j=1}^{d}p_{\lambda_{n,k}}(\beta_{kj}) \right\}.}
#'     In the M step of EM algorithm the
#'     function \deqn{\tilde{Q}(\boldsymbol\Psi,\boldsymbol\Psi^{(m)})
#'     =\sum\limits_{k=1}^{K} \tilde{Q}_{k}(\boldsymbol\Psi_k,
#'     \boldsymbol\Psi^{(m)}_k) =
#'     \sum\limits_{k=1}^{K} \left[{Q}_{k}(\boldsymbol\Psi_k,
#'     \boldsymbol\Psi^{(m)}_k) - \pi_{k}^\alpha\left\{
#'     \sum\limits_{j=1}^{d}p_{\lambda_{n,k}}(\beta_{kj}) \right\}\right]}
#'     is maximized. Since the penalty function is singular at origin, we
#'     use a local quadratic approximation (LQA) for the penalty as
#'     follows, \deqn{\mathbf{p}^\ast_{k,\boldsymbol\lambda_{n}}
#'     (\boldsymbol\beta,\boldsymbol\beta^{(m)})
#'     =(\pi_{k}^{(m)})^{\alpha}\sum\limits_{j=1}^{d}\left\{
#'     p_{\lambda_{n,k}}(\beta_{kj}^{(m)}) + { p^{\prime}_{\lambda_{n,k}}
#'     (\beta_{kj}^{(m)})  \over 2\beta_{kj}^{(m)}}(\beta_{kj}^{2} -
#'     {\beta_{kj}^{(m)}}^{2}) \right\}.} Therefore maximizing \eqn{Q} is
#'     equivalent to maximizing the
#'     function \deqn{ {Q}^\ast(\boldsymbol\Psi,\boldsymbol\Psi^{(m)})
#'     =\sum\limits_{k=1}^{K} {Q}^\ast_{k}(\boldsymbol\Psi_k,
#'     \boldsymbol\Psi^{(m)}_k) = \sum\limits_{k=1}^{K}
#'     \left[{Q}_{k}(\boldsymbol\Psi_k,\boldsymbol\Psi^{(m)}_k)-
#'     \mathbf{p}^\ast_{k,\boldsymbol\lambda_{n}}(\boldsymbol\beta,
#'     \boldsymbol\beta^{(m)})\right].}
#'     In case of Log-Normal sub-distributions, the maximizers of \eqn{Q_k}
#'     functions are as follows. Given the data and current estimates of
#'     parameters, the maximizers are \deqn{{\boldsymbol\beta}^{(m+1)}_{k}
#'     =({\boldsymbol z}^{\prime}\boldsymbol\tau^{(m)}_{k}{\boldsymbol z}+
#'     \varpi_{k}(\boldsymbol\beta_{kj}^{(m)}))^{-1}{\boldsymbol z}^{\prime}
#'     \boldsymbol\tau^{(m)}_{k}T^{(m)}_{k},}
#'     where \eqn{\varpi_{k}(\boldsymbol\beta_{kj}^{(m)})={diag}
#'     \left(\left(\pi_{k}^{(m+1)}\right)^\alpha
#'     \frac{{p}^{\prime}_{\lambda_{n},k}(\boldsymbol\beta_{kj}^{(m)})}
#'     {\boldsymbol\beta_{kj}^{(m)}}\right)}
#'     and \eqn{\sigma_{k}^{(m+1)}} is equal to \deqn{\sigma_{k}^{(m+1)}
#'     =\sqrt{\frac{\sum\limits_{i=1}^{n}\tau^{(m)}_{ik} (t^{(m)}_{ik}
#'     -{\boldsymbol z}_{i}\boldsymbol\beta^{(m)}_{k})^{2}}
#'     {\sum\limits_{i=1}^{n}\tau^{(m)}_{ik} {\left[\delta_{i}
#'     +(1-\delta_{i})\{A(w^{(m)}_{ik})[A(w^{(m)}_{ik})-
#'     w^{(m)}_{ik}]\}\right]}}}.}
#'     For the Weibull distribution, on the other hand,  we have
#'     \eqn{\tilde{\boldsymbol\Psi}^{(m+1)}_k
#'     =\tilde{\boldsymbol\Psi}^{(m)}_k
#'     - 0.5^{\kappa}\left[{H_{k}^{p,(m)}}\right]^{-1}I_{k}^{p,(m)}},
#'     where \eqn{H^p_{k}=H_k+h(\boldsymbol\Psi_k)}
#'     is the penalized version of hessian matrix
#'     and \eqn{I^p_{k}=I_k+h(\boldsymbol\Psi_k)\boldsymbol\Psi_k}
#'     is the penalized version of vector of first derivatives evaluated
#'     at \eqn{\tilde{\boldsymbol\Psi}_k^{(m)}}.
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S.
#' (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{fmrsfit-class}} object which includes parameter
#'     estimates of an FMRs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' n = 500
#' deviance = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrsgendata(n = n, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), deviance = deviance,
#'                      mixProp =mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrsmle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDeviance = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#'
#' res.lam <- fmrstunsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                       nComp = ncomp(res.mle), disFamily = "lnorm",
#'                       initCoeff=c(coefficients(res.mle)),
#'                       initDeviance = deviance(res.mle),
#'                       initmixProp = mixProp(res.mle),
#'                       penFamily = "adplasso")
#' res.var <- fmrsvarsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                       nComp = ncomp(res.mle), disFamily = "lnorm",
#'                       initCoeff=c(coefficients(res.mle)),
#'                       initDeviance = deviance(res.mle),
#'                       initmixProp = mixProp(res.mle),
#'                       penFamily = "adplasso",
#'                       lambPen = slot(res.lam, "lambPen"))
#'
#' coefficients(res.var)[-1,]
#' round(coefficients(res.var)[-1,],5)
#' @export
fmrsvarsel <- function(y,
                       x,
                       delta,
                       nComp,
                       disFamily = "lnorm",
                       initCoeff,
                       initDeviance,
                       initmixProp,
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
  if(is.null(initCoeff) | is.null(initDeviance) | is.null(initmixProp))
    stop("Initial values are not specified.")
  if(length(initCoeff) != nComp*nCov+nComp | length(initmixProp)!=nComp |
     length(initDeviance)!=nComp)
    stop("The length of initial values are not correctly specified.")
  if(!is.matrix(x))
    stop("Provide a matix for covariates.")
  nCov = dim(x)[2]
  n = length(y)
  if(dim(x)[1]!=n)
    stop("The length of observations and rows of design matrix
         does not match.")

  coef0 <- matrix(initCoeff, nrow = nComp, ncol = nCov+1, byrow = TRUE)

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
    model = "FMR"
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
           Initial.Deviance = as.double(initDeviance),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Deviance.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
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
    model = "FMAFTR"
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
           Initial.Deviance = as.double(initDeviance),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Deviance.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
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
    model = "FMAFTR"
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
           Initial.Deviance = as.double(initDeviance),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(convepsNR),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Deviance.Hat = as.double(rep(0,nComp)),
           mixProp.Hat = as.double(rep(0,nComp)),
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

  fit <- new("fmrsfit", y = y,
             delta = delta,
             x = x,
             nobs = n,
             ncov = nCov,
             ncomp = nComp,
             coefficients = array(rbind(res$Intecept.Hat,
                                        matrix(res$Coefficient.Hat,
                                               nrow = nCov, byrow = FALSE)),
                                  dim = c(nCov+1, nComp),
                                  dimnames = list(xnames,comnames)),
             deviance = array(res$Deviance.Hat, dim = c(1,nComp),dimnames =
                                list(NULL,comnames)),
             mixProp = array(res$mixProp.Hat, dim = c(1,nComp),dimnames =
                               list(NULL,comnames)),
             logLik = res$LogLikelihood,
             BIC = res$BIC,
             nIterEMconv = res$Max.iterEM.used,
             disFamily = disFamily,
             penFamily = penFamily,
             lambPen = array(lambPen, dim = c(1,nComp),dimnames =
                               list(NULL,comnames)),
             model = model,
             fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
                            dim = c(n, nComp), dimnames =
                              list(NULL,comnames)),
             residuals = array(matrix(res$residual, nrow =n, byrow = FALSE),
                               dim = c(n, nComp), dimnames =
                                 list(NULL,comnames)),
             weights = array(matrix(res$tau, nrow = n, byrow = FALSE),
                             dim = c(n, nComp), dimnames =
                               list(NULL,comnames))
  )
  return(fit)
}

