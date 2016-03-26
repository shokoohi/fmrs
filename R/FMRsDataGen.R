#' @title Generating data from FMRs models
#'
#' @name fmrs.gen.data
#' @description This function will generate a data set from Finite Mixture of AFT regression models or Finite Mixture of Regression models.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @return A list including a vector of observations \code{y}, a vector of censoring indicators \code{delta} and a matrix of covariates \code{x}
#' @param disFamily Specify the family of sub-distributioons. The options are \code{"lnormal"} for Log-Normal, \code{"norm"} for Normal and \code{"weibull"} for Weibull.
#' @param n A numeric value represents number of observations (sample size)
#' @param nComp A numeric value represents the order (number of components) of an FMRs model
#' @param nCov A numberic value represents the number of covariates in design matrix
#' @param coeff A vector of all regression coefficients including intercepts. It must be a vector of size \code{nComp} by \code{nCov+1}.
#' @param deviance A vector of positive values for dispersion parameters of sub-distributions in FMRs models
#' @param pi A vector of mixing proportions which their sum must be one
#' @param rho A numeric value in [-1, 1] which represents the correlation between covariates of design matrix
#' @param umax A numeric value represents the upper bound in Uniform distribution for censoring
#' @import stats
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S. (2016 submitted) Variable Selection in Mixture of Survival Models
#' @keywords FMR, AFT, FMRs, Normal, Log-Normal, Weibull
#' @concept fmr, aft
#'@examples
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
#' @export
fmrs.gen.data <- function(n,
                          nComp,
                          nCov,
                          coeff,
                          deviance,
                          pi,
                          rho,
                          umax,
                          disFamily = "lnorm"
)
{
  if(sum(pi) != 1)
    stop("The sum of mixing proportions must be 1.")
  if(sum(deviance <= 0) != 0)
    stop("Dispersion parameters cannot be zero or negative.")
  if(rho > 1 | rho < -1)
    stop("The correlation cannot be less than -1 or greater thatn 1.")

  mu <- rep(0, nCov)
  Sigma <- diag(nCov)

  for(i in 1:nCov){
    for(j in 1:nCov){
      Sigma[i,j] <- rho^abs(i-j)
    }}

  X <- matrix(rnorm(nCov * n), n)
  X <- scale(X, TRUE, FALSE)
  X <- X %*% svd(X, nu = 0)$v
  X <- scale(X, FALSE, TRUE)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), nCov) %*% t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    cX = drop(X)
  else cX = t(X)
  cX <- scale(cX)
  colnames(cX) <-  paste("X", 1:nCov,sep = ".")

  coef0 <- matrix(coeff, nrow = nComp, ncol = nCov+1, byrow = T)
  pi0 <- cumsum(pi)

  yobs <-c()
  c <- rep()
  dlt <- c()
  u <- c()
  tobs <- c()

  if(disFamily == "lnorm"){
    for(i in 1:n){
      epss <- rnorm(1)
      u1 <- runif(1)
      k = length(which(pi0<=u1)) + 1
      u[i] = k
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + deviance[k] * epss

      c[i] <- log(runif(1, 0, umax))
      tobs[i] <- exp(min(yobs[i],c[i]))
      dlt[i] <- (yobs[i] < c[i])*1
    }
  }else if(disFamily=="norm"){
    for(i in 1:n){
      epss <- rnorm(1)
      u1 <- runif(1)
      k = length(which(pi0<=u1)) + 1
      u[i] = k
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + deviance[k] * epss
      tobs[i] <- yobs[i]
      dlt[i] <- 1
    }
  }else if(disFamily=="weibull"){
    for(i in 1:n){
      ext <- log(rexp(1))
      u1 <- runif(1)
      k = length(which(pi0<=u1)) + 1
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + deviance[k] * ext

      c[i]<- log(runif(1, 0, umax))
      tobs[i] <- exp(min(yobs[i],c[i]))
      dlt[i] <- (yobs[i] < c[i])*1
    }
  } else{
    stop("The family of sub-distributions are not specified correctly.")
  }
  return(list(y = tobs, delta = dlt, x = cX, disFamily = disFamily))
}

