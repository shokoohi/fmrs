#' @rdname weights-methods
#' @rdname fmrsfit-class
#' @aliases weights,weights-method
setMethod("weights", signature = "fmrsfit", weights.fmrsfit)

#' @rdname residuals-methods
#' @rdname fmrsfit-class
#' @aliases residuals,residuals-method
setMethod("residuals", signature = "fmrsfit", residuals.fmrsfit)

#' @rdname nobs-methods
#' @rdname fmrsfit-class
#' @aliases nobs,nobs-method
setMethod("nobs", signature = "fmrsfit", nobs.fmrsfit)

#' @rdname ncov-methods
#' @rdname fmrsfit-class
#' @aliases ncov,ncov-method
setMethod("ncov", signature = "fmrsfit", ncov.fmrsfit)

#' @rdname ncomp-methods
#' @rdname fmrsfit-class
#' @aliases ncomp,ncomp-method
setMethod("ncomp", signature = "fmrsfit", ncomp.fmrsfit)

#' @rdname mixProp-methods
#' @rdname fmrsfit-class
#' @aliases mixProp,mixProp-method
setMethod("mixProp", signature = "fmrsfit", mixProp.fmrsfit)

#' @rdname logLik-methods
#' @rdname fmrsfit-class
#' @aliases logLik,logLik-method
setMethod("logLik", signature = "fmrsfit", logLik.fmrsfit)

#' @rdname fitted-methods
#' @rdname fmrsfit-class
#' @aliases fitted,fitted-method
setMethod("fitted", signature = "fmrsfit", fitted.fmrsfit)

#' @rdname dispersion-methods
#' @rdname fmrsfit-class
#' @aliases dispersion,dispersion-method
setMethod("dispersion", signature = "fmrsfit", dispersion.fmrsfit)

#' @rdname coefficients-methods
#' @rdname fmrsfit-class
#' @aliases coefficients,coefficients-method
setMethod("coefficients", signature = "fmrsfit", coefficients.fmrsfit)

#' @rdname BIC-methods
#' @rdname fmrsfit-class
#' @aliases BIC,BIC-method
setMethod("BIC", signature = "fmrsfit", BIC.fmrsfit)

#' @rdname summary-methods
#' @rdname fmrsfit-class
#' @aliases summary,summary-method
setMethod("summary", signature = "fmrsfit", summary.fmrsfit)

#' @rdname summary-methods
#' @rdname fmrstunpar-class
#' @aliases summary,summary-method
setMethod("summary", signature = "fmrstunpar", summary.fmrstunpar)

#' @rdname show-methods
#' @rdname fmrsfit-class
#' @aliases show,show-method
setMethod("show", signature = "fmrsfit", show.fmrsfit)

#' @rdname show-methods
#' @rdname fmrstunpar-class
#' @aliases show,show-method
setMethod("show", signature = "fmrstunpar", show.fmrstunpar)

#' @rdname fmrs.mle-methods
#' @aliases fmrs.mle-method
setMethod(f="fmrs.mle", definition=function(y,
                                         delta,
                                         x,
                                         nComp = 2,
                                         disFamily = "lnorm",
                                         initCoeff,
                                         initDispersion,
                                         initmixProp,
                                         lambRidge = 0,
                                         nIterEM = 400,
                                         nIterNR = 2,
                                         conveps = 1e-8,
                                         convepsEM = 1e-8,
                                         convepsNR = 1e-8,
                                         porNR = 2,
                                         activeset){
  if(missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
  if(missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
  if(missing(delta) & (disFamily!="norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

    if((nComp<2) ) {
    stop("An interger greater than 1 for the order of mixture model
         must be provided.")
  }
  nCov = dim(x)[2]
  n = length(y)
  if(missing(initCoeff)) initCoeff = rnorm((nCov+1)*nComp)
  if(missing(initDispersion)) initDispersion = rep(1,nComp)
  if(missing(initmixProp)) initmixProp = rep(1/nComp,nComp)

  if(missing(activeset)) activeset = matrix(1, nrow = nCov+1, ncol = nComp)
  if(any(activeset!=0 & activeset!=1))
    stop("activeset must be a matrix of dimention nCov+1 by nComp with
         only 0 and 1 values.")

    coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("X",1:nCov,sep=".")))
  } else{
    xnames <- c("Intercept",colnames(x))
  }
  comnames <- c(paste("Comp",1:nComp,sep="."))

  if(disFamily == "norm"){
    model = "FMR"
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
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
           tau = as.double(rep(0,n*nComp)),
           actset = as.integer(activeset),
           disnorm = as.integer(1)
    )
  }else if(disFamily == "lnorm"){
    model = "FMAFTR"
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
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
           tau = as.double(rep(0,n*nComp)),
           actset = as.integer(activeset),
           disnorm = as.integer(2)
    )

  }else if(disFamily == "weibull"){
    model = "FMAFTR"
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps.em = as.double(convepsEM),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
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
           tau = as.double(rep(0,n*nComp)),
           actset = as.integer(activeset)
    )

    if(res$LogLikelihood=='NaN'){
      fmrs.mle2 <- function(y,
                            x,
                            delta,
                            Lambda.Ridge,
                            Num.Comp,
                            Num.Cov,
                            Sample.Size,
                            Num.iterationEM,
                            Num.iterationNR,
                            PortionNF,
                            Initial.Intercept,
                            Initial.Coefficient,
                            Initial.Dispersion,
                            Initial.mixProp,
                            conv.eps.em,
                            actset
      )
      {
        require(survival)
        x <- matrix(c(x),Sample.Size,Num.Cov)
        acs=matrix(c(actset),nCov+1,nComp)
        new_pi0 <- pi0 <- c(Initial.mixProp)
        new_beta0 <- beta0 <- matrix(Initial.Coefficient,nCov,nComp)
        new_alpha0 <- alpha0 <- c(Initial.Intercept)
        new_sigma0 <- sigma0 <- c(Initial.Dispersion)

        phi <- matrix(0, Sample.Size, Num.Comp)
        W <- predict <- residual <- matrix(0, Sample.Size, Num.Comp)
        sumwi <- matrix(0, Num.Comp)

        emiter = 0
        CONV = 0
        while ((emiter<Num.iterationEM) & (CONV!=1)) {
          for (k1 in 1:Num.Comp) {
            sumwi[k1] = 0
          }

          for (i in 1:Sample.Size) {
            sumi = 0.0
            for (k1 in 1:Num.Comp) {
              mui = sum(x[i,] * beta0[,k1] * acs[-1,k1])
              mui = mui + alpha0[k1] * acs[1,k1]
              deni = (((1 / sigma0[k1])* exp((log(y[i]) - mui) / sigma0[k1]))^(delta[i])) * exp(-exp((log(y[i]) - mui) / sigma0[k1]))
              if(deni==0)
                deni= 0.000001
              phi[i,k1] = pi0[k1] * deni
              sumi = sumi + phi[i,k1]
            }
            for(k1 in 1:Num.Comp)
            {
              W[i,k1] = phi[i,k1] / sumi
              if(W[i,k1]<1e-10) W[i,k1] = 1e-10
              sumwi[k1] = sumwi[k1] + W[i,k1]
            }
          }

          for(k1 in 1:Num.Comp){
            new_pi0[k1] = sumwi[k1] / Sample.Size;
          }

          for(k1 in 1:Num.Comp)
          {
            newX = x[,acs[-1,k1]==1]
            if(acs[1,k1]==1){
              res <- survival::survreg(Surv(time = y, event = delta, type = c('right')) ~ 1 + newX , weights = W[,k1],
                             dist="weibull", init=c(alpha0[k1],beta0[acs[-1,k1]==1,k1]), scale=sigma0[k1],
                             control=list(maxiter=Num.iterationNR, rel.tolerance=conv.eps.em,
                                                     toler.chol=conv.eps.em, iter.max=Num.iterationNR, debug=0, outer.max=10),parms=NULL,model=FALSE, x=FALSE,
                             y=TRUE, robust=FALSE, score=FALSE)
              new_alpha0[k1] <- as.double(coef(res)[1])
              new_beta0[acs[-1,k1]==1,k1] <- as.double(coef(res)[-1])
              new_beta0[acs[-1,k1]==0,k1] <- 0
            }
            else{
              res <- survival::survreg(Surv(time = y, event = delta, type = c('right')) ~ -1 + newX , weights = W[,k1],
                             dist="weibull", init=c(beta0[acs[-1,k1]==1,k1]), scale=sigma0[k1],
                             control=list(maxiter=Num.iterationNR, rel.tolerance=conv.eps.em,
                                                     toler.chol=conv.eps.em, iter.max=Num.iterationNR, debug=0, outer.max=10),parms=NULL,model=FALSE, x=FALSE,
                             y=TRUE, robust=FALSE, score=FALSE)
              new_alpha0[k1] <- 0
              new_beta0[acs[-1,k1]==1,k1] <- as.double(coef(res))
              new_beta0[acs[-1,k1]==0,k1] <- 0
            }
            sumi3 = 0.0;
            sumi5 = 0.0;

            for(i in 1:Sample.Size){
              mui = 0.0;
              mui = sum(x[i,] * new_beta0[,k1] * acs[-1,k1])
              mui = mui + new_alpha0[k1] * acs[1,k1]
              sumi3 = sumi3 + W[i,k1] * (- delta[i] / sigma0[k1] + ((log(y[i]) - mui) / (sigma0[k1] * sigma0[k1])) * ( exp( (log(y[i]) - mui) / sigma0[k1]) - delta[i])   )
              sumi5 = sumi5 + W[i,k1] * ( delta[i] / (sigma0[k1] * sigma0[k1]) +  ((log(y[i]) - mui) / (sigma0[k1] * sigma0[k1] * sigma0[k1])) * (2 * delta[i] - (2 + (log(y[i]) - mui) / sigma0[k1]) * exp( (log(y[i]) - mui) / sigma0[k1])  ))
            }
            new_sigma0[k1] = sigma0[k1] - (1 / sumi5) * sumi3;
          }

          emiter = emiter + 1

          diff = sum((new_sigma0 - sigma0)^2) + sum((new_alpha0 - alpha0)^2) + sum((new_beta0 - beta0)^2) + sum((new_pi0-pi0)^2)
          if(diff<conv.eps.em)
            CONV = 1
          sigma0 = new_sigma0
          pi0 = new_pi0
          alpha0 = new_alpha0
          beta0 = new_beta0
        }

        loglike1 = 0.0


        for(i in 1:Sample.Size){
          sumi = 0.0;
          for(k1 in 1:Num.Comp){
            mui = sum(x[i,] * new_beta0[,k1] * acs[-1,k1]) + new_alpha0[k1] * acs[1,k1]
            deni = (((1 / new_sigma0[k1])* exp((log(y[i]) - mui) / new_sigma0[k1]))^(delta[i])) * exp(-exp((log(y[i]) - mui) / new_sigma0[k1]))
            phi[i,k1] = new_pi0[k1] * deni
            sumi = sumi + phi[i,k1]
          }
          loglike1 = loglike1 + log(sumi)
        }

        BIC = loglike1 - 0.5 * Num.Comp * Num.Cov * log(Sample.Size)
        EBIC5 = loglike1 - 0.5 * Num.Comp * Num.Cov * log(Sample.Size) - 0.5 * (Num.Comp * Num.Cov) * log(Num.Cov)
        EBIC1 = loglike1 - 0.5 * (Num.Comp * Num.Cov) * log(Sample.Size) - (Num.Comp * Num.Cov) * log(Num.Cov)
        AIC = loglike1 - (Num.Comp * Num.Cov)
        GCV = (loglike1) / (Sample.Size * (1 - Num.Comp * Num.Cov / Sample.Size)^ 2)
        GIC = loglike1 - 0.5 * (Num.Comp * Num.Cov) * log(Sample.Size)
        MaxEMiter = emiter

        for (i in 1:Sample.Size) {
          sumi = 0.0
          for (k1 in 1:Num.Comp) {
            mui = sum(x[i,] * new_beta0[,k1] * acs[-1,k1]) + new_alpha0[k1] * acs[1,k1]
            deni = (((1 / new_sigma0[k1])* exp((log(y[i]) - mui) / new_sigma0[k1]))^(delta[i])) * exp(-exp((log(y[i]) - mui) / new_sigma0[k1]))
            if(deni==0)
              deni= 0.000001
            phi[i,k1] = new_pi0[k1] * deni
            sumi = sumi + phi[i,k1]
          }
          for(k1 in 1:Num.Comp)
          {
            W[i,k1] = phi[i,k1] / sumi
          }
        }


        for(k1 in 1:Num.Comp){
          for(i in 1:Sample.Size){
            mui = sum(x[i,] * new_beta0[,k1] * acs[-1,k1]) + new_alpha0[k1] * acs[1,k1]
            predict[i,k1] = exp(mui)
            residual[i,k1] = y[i] - exp(mui)
          }
        }


        list(Intecept.Hat = c(new_alpha0),
             Coefficient.Hat = c(new_beta0),
             Dispersion.Hat = c(new_sigma0),
             mixProp.Hat = c(new_pi0),
             LogLikelihood = loglike1,
             BIC = BIC,
             AIC = AIC,
             GCV = GCV,
             EBIC1 = EBIC1,
             EBIC5 = EBIC5,
             GIC = GIC,
             predict = c(predict),
             residual = c(residual),
             tau = c(W),
             Max.iterEM.used = MaxEMiter
        )
      }

      res <- fmrs.mle2(
             y = as.double(y),
             x = c(as.double(as.vector(unlist(x)))),
             delta = as.double(delta),
             Lambda.Ridge = as.double(lambRidge),
             Num.Comp = as.integer(nComp),
             Num.Cov = as.integer(nCov),
             Sample.Size = as.integer(n),
             Num.iterationEM = as.integer(nIterEM),
             Num.iterationNR = 1,
             PortionNF = as.integer(porNR),
             Initial.Intercept = as.double(c(coef0[,1])),
             Initial.Coefficient = as.double(c(t(coef0[,-1]))),
             Initial.Dispersion = as.double(initDispersion),
             Initial.mixProp = as.double(initmixProp),
             conv.eps.em = as.double(convepsEM),
             actset = as.integer(activeset)
      )


      print("The results are based on using survreg in survival package.")
    }

  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }


  fit <- new("fmrsfit",
             y = y,
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
             dispersion = array(res$Dispersion.Hat, dim =
                                c(1,nComp),dimnames = list(NULL,comnames)),
             mixProp = array(res$mixProp.Hat, dim =
                               c(1,nComp),dimnames = list(NULL,comnames)),
             logLik = res$LogLikelihood,
             BIC = res$BIC,
             nIterEMconv = res$Max.iterEM.used,
             disFamily = disFamily,
             lambRidge = lambRidge,
             model = model,
             fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
                            dim = c(n, nComp), dimnames =
                              list(NULL,comnames)),
             residuals = array(matrix(res$residual, nrow =n, byrow = FALSE),
                               dim = c(n, nComp),
                               dimnames = list(NULL,comnames)),
             weights = array(matrix(res$tau, nrow = n, byrow = FALSE),
                             dim = c(n, nComp), dimnames =
                               list(NULL,comnames)),
             activeset = array(matrix(c(activeset), nrow=nCov+1, byrow = FALSE),
                              dim = c(nCov+1, nComp),
                              dimnames = list(xnames,comnames))
  )
  return(fit)
})




#' @rdname fmrs.tunsel-methods
#' @aliases fmrs.tunsel-method
setMethod(f="fmrs.tunsel", definition=function(y,
                                               delta,
                                               x,
                                               nComp,
                                               disFamily = "lnorm",
                                               initCoeff,
                                               initDispersion,
                                               initmixProp,
                                               penFamily = "lasso",
                                               lambRidge = 0,
                                               nIterEM = 2000,
                                               nIterNR = 2,
                                               conveps = 1e-8,
                                               convepsEM = 1e-8,
                                               convepsNR = 1e-8,
                                               porNR = 2,
                                               gamMixPor = 1,
                                               activeset,
                                               lambMCP,
                                               lambSICA
                                               ){
  if(missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
  if(missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
  if(missing(delta) & (disFamily!="norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

  if(missing(lambMCP))
    lambMCP = 2/(1-max(cor(x)[cor(x)!=1]))

  if(missing(lambSICA))
    lambSICA = 5.0

  if((nComp<2) ) {
    stop("An interger greater than 2 for the order of mixture model
         must be provided.")
  }
  nCov = dim(x)[2]
  n = length(y)
  if(missing(initCoeff)) initCoeff = rnorm((nCov+1)*nComp)
  if(missing(initDispersion)) initDispersion = rep(1,nComp)
  if(missing(initmixProp)) initmixProp = rep(1/nComp,nComp)
  if(missing(activeset)) activeset = matrix(1, nrow = nCov+1, ncol = nComp)
  if(any(activeset!=0 & activeset!=1))
    stop("activeset must be a matrix of dimention nCov+1 by nComp with
         only 0 and 1 values.")

  coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("X",1:nCov,sep=".")))
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
  else {
    stop("Penalty is not correctly specified.")
    }

  if(disFamily == "norm"){
    model = "FMR"
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp)),
           actset = as.integer(activeset),
           tuneGam1 = as.double(lambMCP),
           tuneGam1 = as.double(lambSICA)
    )
  }else if(disFamily == "lnorm"){
    model = "FMAFTR"
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp)),
           actset = as.integer(activeset),
           tuneGam1 = as.double(lambMCP),
           tuneGam1 = as.double(lambSICA)
    )
  }else if(disFamily == "weibull"){
    model = "FMAFTR"
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           Num.NRiteration = as.double(nIterNR),
           Num.PortionNF = as.double(porNR),
           conv.eps = as.double(convepsNR),
           GamMixPortion = as.double(gamMixPor),
           Opt.Lambda = as.double(rep(0,nComp)),
           actset = as.integer(activeset),
           tuneGam1 = as.double(lambMCP),
           tuneGam1 = as.double(lambSICA)
    )
  }else{
    stop("The family of sub-distributions is not specified correctly.")
  }


  lambdafit <- new("fmrstunpar",
                   ncomp = nComp,
                   lambPen = array(res$Opt.Lambda, dim = c(1,nComp),
                                   dimnames = c(list(NULL,
                                                     c(paste("Comp", 1:nComp,
                                                             sep = "."))))),
                   lambRidge = lambRidge,
                   disFamily = disFamily,
                   penFamily = penFamily,
                   MCPGam = lambMCP,
                   SICAGam = lambSICA,
                   model = model,
                   activeset = array(matrix(c(activeset), nrow=nCov+1, byrow = FALSE),
                                     dim = c(nCov+1, nComp),
                                     dimnames = list(xnames,comnames))
  )
  return(lambdafit)
}
)

#' @rdname fmrs.varsel-methods
#' @aliases fmrs.varsel-method
setMethod(f="fmrs.varsel", definition=function(y,
                                               delta,
                                               x,
                                               nComp,
                                               disFamily = "lnorm",
                                               initCoeff,
                                               initDispersion,
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
                                               gamMixPor = 1,
                                               activeset,
                                               lambMCP,
                                               lambSICA
                                               ){
  if(missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
  if(missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
  if(missing(delta) & (disFamily!="norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

  if(missing(lambMCP))
    lambMCP = 2/(1-max(cor(x)[cor(x)!=1]))

  if(missing(lambSICA))
    lambSICA = 5.0

  if((nComp<2) ) {
    stop("An interger greater than 2 for the order of mixture model
         must be provided.")
  }
  nCov = dim(x)[2]
  n = length(y)
  if(missing(initCoeff)) initCoeff = rnorm((nCov+1)*nComp)
  if(missing(initDispersion)) initDispersion = rep(1,nComp)
  if(missing(initmixProp)) initmixProp = rep(1/nComp,nComp)

  if(missing(activeset)) activeset = matrix(1, nrow = nCov+1, ncol = nComp)
  if(any(activeset!=0 & activeset!=1))
    stop("activeset must be a matrix of dimention nCov+1 by nComp with
         only 0 and 1 values.")

  coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

  if(is.null(colnames(x))){
    xnames <- c("Intercept",c(paste("X",1:nCov,sep=".")))
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
  else {stop("Penalty is not correctly specified.")  }

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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
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
           tau = as.double(rep(0,n*nComp)),
           actset = as.integer(activeset),
           tuneGam1 = as.double(lambMCP),
           tuneGam1 = as.double(lambSICA)
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(conveps),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
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
           tau = as.double(rep(0,n*nComp)),
           actset = as.integer(activeset),
           tuneGam1 = as.double(lambMCP),
           tuneGam1 = as.double(lambSICA)
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
           Initial.Dispersion = as.double(initDispersion),
           Initial.mixProp = as.double(initmixProp),
           conv.eps = as.double(convepsNR),
           conv.eps.em = as.double(convepsEM),
           GamMixPortion = as.double(gamMixPor),
           Intecept.Hat = as.double(rep(0,nComp)),
           Coefficient.Hat = as.double(rep(0,nComp*nCov)),
           Dispersion.Hat = as.double(rep(0,nComp)),
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
           tau = as.double(rep(0,n*nComp)),
           actset = as.integer(activeset),
           tuneGam1 = as.double(lambMCP),
           tuneGam1 = as.double(lambSICA)
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
             dispersion = array(res$Dispersion.Hat, dim = c(1,nComp),dimnames =
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
             MCPGam = lambMCP,
             SICAGam = lambSICA,
             model = model,
             fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
                            dim = c(n, nComp), dimnames =
                              list(NULL,comnames)),
             residuals = array(matrix(res$residual, nrow =n, byrow = FALSE),
                               dim = c(n, nComp), dimnames =
                                 list(NULL,comnames)),
             weights = array(matrix(res$tau, nrow = n, byrow = FALSE),
                             dim = c(n, nComp), dimnames =
                               list(NULL,comnames)),
             activeset = array(matrix(c(activeset), nrow=nCov+1, byrow = FALSE),
                               dim = c(nCov+1, nComp),
                               dimnames = list(xnames,comnames))
  )
  return(fit)
}

)


#' @rdname fmrs.gendata-methods
#' @aliases fmrs.gendata-method
setMethod(f="fmrs.gendata", definition=function(nObs,
                                                nComp,
                                                nCov,
                                                coeff,
                                                dispersion,
                                                mixProp,
                                                rho,
                                                umax,
                                                disFamily = "lnorm")

{
  if(missing(disFamily)) disFamily = "lnorm"

  if(sum(mixProp) != 1)
    stop("The sum of mixing proportions must be 1.")
  if(sum(dispersion <= 0) != 0)
    stop("Dispersion parameters cannot be zero or negative.")
  if(rho > 1 | rho < -1)
    stop("The correlation cannot be less than -1 or greater thatn 1.")

  mu <- rep(0, nCov)
  Sigma <- diag(nCov)

  for(i in 1:nCov){
    for(j in 1:nCov){
      Sigma[i,j] <- rho^abs(i-j)
    }}

  X <- matrix(rnorm(nCov * nObs), nObs)
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
  if (nObs == 1)
    cX = drop(X)
  else cX = t(X)
  cX <- scale(cX)
  colnames(cX) <-  paste("X", 1:nCov,sep = ".")

  coef0 <- matrix(coeff, nrow = nComp, ncol = nCov+1, byrow = TRUE)
  mixProp0 <- cumsum(mixProp)

  yobs <-c()
  c <- rep()
  dlt <- c()
  u <- c()
  tobs <- c()

  if(disFamily == "lnorm"){
    for(i in 1:nObs){
      epss <- rnorm(1)
      u1 <- runif(1)
      k = length(which(mixProp0<=u1)) + 1
      u[i] = k
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * epss

      c[i] <- log(runif(1, 0, umax))
      tobs[i] <- exp(min(yobs[i],c[i]))
      dlt[i] <- (yobs[i] < c[i])*1
    }
  }else if(disFamily=="norm"){
    for(i in 1:nObs){
      epss <- rnorm(1)
      u1 <- runif(1)
      k = length(which(mixProp0<=u1)) + 1
      u[i] = k
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * epss
      tobs[i] <- yobs[i]
      dlt[i] <- 1
    }
  }else if(disFamily=="weibull"){
    for(i in 1:nObs){
      ext <- log(rexp(1))
      u1 <- runif(1)
      k = length(which(mixProp0<=u1)) + 1
      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * ext

      c[i]<- log(runif(1, 0, umax))
      tobs[i] <- exp(min(yobs[i],c[i]))
      dlt[i] <- (yobs[i] < c[i])*1
    }
  } else{
    stop("The family of sub-distributions are not specified correctly.")
  }
  return(list(y = tobs, delta = dlt, x = cX, disFamily = disFamily))
}
)
