#' @rdname weights-methods
#' @aliases weights,weights-method
setMethod("weights", signature = "fmrsfit", weights.fmrsfit)

#' @rdname residuals-methods
#' @aliases residuals,residuals-method
setMethod("residuals", signature = "fmrsfit", residuals.fmrsfit)

#' @rdname nobs-methods
#' @aliases nobs,nobs-method
setMethod("nobs", signature = "fmrsfit", nobs.fmrsfit)

#' @rdname ncov-methods
#' @aliases ncov,ncov-method
setMethod("ncov", signature = "fmrsfit", ncov.fmrsfit)

#' @rdname ncomp-methods
#' @aliases ncomp,ncomp-method
setMethod("ncomp", signature = "fmrsfit", ncomp.fmrsfit)

#' @rdname mixProp-methods
#' @aliases mixProp,mixProp-method
setMethod("mixProp", signature = "fmrsfit", mixProp.fmrsfit)

#' @rdname logLik-methods
#' @aliases logLik,logLik-method
setMethod("logLik", signature = "fmrsfit", logLik.fmrsfit)

#' @rdname fitted-methods
#' @aliases fitted,fitted-method
setMethod("fitted", signature = "fmrsfit", fitted.fmrsfit)

#' @rdname dispersion-methods
#' @aliases dispersion,dispersion-method
setMethod("dispersion", signature = "fmrsfit", dispersion.fmrsfit)

#' @rdname coefficients-methods
#' @aliases coefficients,coefficients-method
setMethod("coefficients", signature = "fmrsfit", coefficients.fmrsfit)

#' @rdname BIC-methods
#' @aliases BIC,BIC-method
setMethod("BIC", signature = "fmrsfit", BIC.fmrsfit)

#' @rdname summary-methods
#' @aliases summary,summary-method
setMethod("summary", signature = "fmrsfit", summary.fmrsfit)

#' @rdname summary-methods
#' @aliases summary,summary-method
setMethod("summary", signature = "fmrstunpar", summary.fmrstunpar)


#' @rdname fmrs.mle-methods
#' @aliases fmrs.mle-method
setMethod(f = "fmrs.mle",
    definition = function(y, delta, x, nComp = 2, disFamily = "lnorm",
    initCoeff, initDispersion, initmixProp,
    lambRidge = 0, nIterEM = 400, nIterNR = 2,
    conveps = 1e-08, convepsEM = 1e-08,
    convepsNR = 1e-08, NRpor = 2, activeset) {
    if (missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
    if (missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
    if (missing(delta) & (disFamily != "norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

    if ((nComp < 2)) {
    stop("An interger greater than 1 for the order of mixture model must be
    provided.")
    }
    if(nComp > 10){
    warning("A very high value is selected for nComp!")
    }

    n = length(y)
    nCov = dim(x)[2]
    if (missing(initCoeff))
    initCoeff = rnorm((nCov + 1) * nComp)
    if (missing(initDispersion))
    initDispersion = rep(1, nComp)
    if (missing(initmixProp))
    initmixProp = rep(1/nComp, nComp)

    if (missing(activeset))
    activeset = matrix(1, nrow = nCov + 1, ncol = nComp)
    if (any(activeset != 0 & activeset != 1))
    stop("activeset must be a matrix of dimention nCov+1 by nComp with only 0
    and 1 values.")

    coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov + 1, byrow = TRUE)

    if (is.null(colnames(x))) {
    xnames <- c("Intercept", c(paste("X", seq_len(nCov), sep = ".")))
    } else {
    xnames <- c("Intercept", colnames(x))
    }
    comnames <- c(paste("Comp", seq_len(nComp), sep = "."))

    if (disFamily == "norm")
    {
    model = "FMR"
    delta = rep(1, n)
    res = .Call(FMR_Norm_MLE,
    as.double(y),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(unlist(coef0[, 1])),
    as.double(unlist(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.double(convepsEM),
    as.integer(1))

    } else if (disFamily == "lnorm")
    {
    model = "FMAFTR"
    logy = log(y)
    res = .Call(FMR_Norm_MLE,
    as.double(logy),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(unlist(coef0[, 1])),
    as.double(unlist(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.double(convepsEM),
    as.integer(2))

    } else if (disFamily == "weibull")
    {
    model = "FMAFTR"
    logy = log(y)
    res = .Call(FMR_Weibl_MLE,
    as.double(logy),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(unlist(coef0[, 1])),
    as.double(unlist(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.double(convepsEM),
    as.integer(nIterNR),
    as.integer(NRpor)
    )

    if (res$LogLik == "NaN") {
    fmrs.mle2 <- function(y, x, delta, Lambda.Ridge, Num.Comp,
    Num.Cov,
    Sample.Size, Num.iterationEM,
    Num.iterationNR,
    PortionNF, Initial.Intercept,
    Initial.Coefficient,
    Initial.Dispersion, Initial.mixProp,
    conv.eps.em, actset) {
    x <- matrix(c(x), Sample.Size, Num.Cov)
    acs = matrix(c(actset), nCov + 1, nComp)
    new_pi0 <- pi0 <- c(Initial.mixProp)
    new_beta0 <- beta0 <- matrix(Initial.Coefficient,
    nCov, nComp)
    new_alpha0 <- alpha0 <- c(Initial.Intercept)
    new_sigma0 <- sigma0 <- c(Initial.Dispersion)

    phi <- matrix(0, Sample.Size, Num.Comp)
    W <- predict <- residual <- matrix(0, Sample.Size,
    Num.Comp)
    sumwi <- matrix(0, Num.Comp)

    emiter = 0
    CONV = 0
    while ((emiter < Num.iterationEM) & (CONV != 1)) {
    sumwi[seq_len(Num.Comp)] = 0
    for (i in seq_len(Sample.Size)) {
    sumi = 0
    for (k1 in seq_len(Num.Comp)) {
    mui = sum(x[i,] * beta0[, k1] *
    acs[-1, k1])
    mui = mui + alpha0[k1] * acs[1, k1]
    deni = (((1/sigma0[k1]) *
    exp((log(y[i]) - mui)/
    sigma0[k1]))^(delta[i])) *
    exp(-exp((log(y[i]) - mui)/sigma0[k1]))
    if (deni == 0)
    deni = 1e-06
    phi[i, k1] = pi0[k1] * deni
    sumi = sumi + phi[i, k1]
    }
    for (k1 in seq_len(Num.Comp)) {
    W[i, k1] = phi[i, k1]/sumi
    if (W[i, k1] < 1e-10)
    W[i, k1] = 1e-10
    sumwi[k1] = sumwi[k1] + W[i, k1]
    }
    }

    for (k1 in seq_len(Num.Comp)) {
    new_pi0[k1] = sumwi[k1]/Sample.Size
    }

    for (k1 in seq_len(Num.Comp)) {
    newX = x[, acs[-1, k1] == 1]
    if (acs[1, k1] == 1) {
    res <- survival::survreg(survival::Surv(time = y, event = delta,
    type = c("right")) ~ 1 + newX,
    weights = W[, k1],
    dist = "weibull",
    init = c(alpha0[k1], beta0[acs[-1, k1] == 1, k1]),
    scale = sigma0[k1],
    control = list(maxiter = Num.iterationNR,
    rel.tolerance = conv.eps.em,
    toler.chol = conv.eps.em,
    iter.max = Num.iterationNR, debug = 0,
    outer.max = 10),
    parms = NULL, model = FALSE, x = FALSE, y = TRUE,
    robust = FALSE, score = FALSE
    )
    new_alpha0[k1] <- as.double(coef(res)[1])
    new_beta0[acs[-1, k1] == 1, k1] <- as.double(coef(res)[-1])
    new_beta0[acs[-1, k1] == 0, k1] <- 0
    } else {
    res <- survival::survreg(survival::Surv(time = y,
    event = delta, type = c("right")) ~ -1 + newX,
    weights = W[, k1], dist = "weibull",
    init = c(beta0[acs[-1, k1] == 1, k1]),
    scale = sigma0[k1],
    control = list(maxiter = Num.iterationNR,
    rel.tolerance = conv.eps.em,
    toler.chol = conv.eps.em, iter.max = Num.iterationNR,
    debug = 0, outer.max = 10), parms = NULL, model = FALSE,
    x = FALSE, y = TRUE, robust = FALSE, score = FALSE)
    new_alpha0[k1] <- 0
    new_beta0[acs[-1, k1] == 1, k1] <- as.double(coef(res))
    new_beta0[acs[-1, k1] == 0, k1] <- 0
    }
    sumi3 = 0
    sumi5 = 0

    for (i in seq_len(Sample.Size)) {
    mui = 0
    mui = sum(x[i,] * new_beta0[, k1] * acs[-1, k1])
    mui = mui + new_alpha0[k1] * acs[1, k1]
    sumi3 = sumi3 + W[i, k1] *
    (-delta[i]/sigma0[k1] + ((log(y[i]) - mui) /
    (sigma0[k1] * sigma0[k1])) *
    (exp((log(y[i]) - mui)/sigma0[k1]) - delta[i]))
    sumi5 = sumi5 + W[i, k1] *
    (delta[i]/(sigma0[k1] * sigma0[k1]) +
    ((log(y[i]) - mui)/
    (sigma0[k1] * sigma0[k1] * sigma0[k1])) *
    (2 * delta[i]-(2 + (log(y[i])-mui)/sigma0[k1]) *
    exp((log(y[i]) - mui)/sigma0[k1])))
    }
    new_sigma0[k1] = sigma0[k1] - (sumi3/sumi5)
    }

    emiter = emiter + 1

    diff = sum((new_sigma0 - sigma0)^2) + sum((new_alpha0 -
    alpha0)^2) +
    sum((new_beta0 - beta0)^2) + sum((new_pi0 - pi0)^2)
    if (diff < conv.eps.em)
    CONV = 1
    sigma0 = new_sigma0
    pi0 = new_pi0
    alpha0 = new_alpha0
    beta0 = new_beta0
    }

    loglike1 = 0

    for (i in seq_len(Sample.Size)) {
    sumi = 0
    for (k1 in seq_len(Num.Comp)) {
    mui = sum(x[i,] * new_beta0[, k1] * acs[-1, k1]) +
    new_alpha0[k1] * acs[1, k1]
    deni = (((1/new_sigma0[k1]) * exp((log(y[i]) - mui) /
    new_sigma0[k1]))^(delta[i])) *
    exp(-exp((log(y[i]) - mui)/new_sigma0[k1]))
    phi[i, k1] = new_pi0[k1] * deni
    sumi = sumi + phi[i, k1]
    }
    loglike1 = loglike1 + log(sumi)
    }

    SumPar = sum(as.integer(activeset))

    BIC = -2 * loglike1 + SumPar * log(Sample.Size)
    AIC = -2*loglike1 + 2*SumPar
    MaxEMiter = emiter

    for (i in seq_len(Sample.Size)) {
    sumi = 0
    for (k1 in seq_len(Num.Comp)) {
    mui = sum(x[i,] * new_beta0[, k1] * acs[-1, k1]) +
    new_alpha0[k1] * acs[1, k1]
    deni = (((1/new_sigma0[k1]) * exp((log(y[i]) - mui) /
    new_sigma0[k1]))^(delta[i])) *
    exp(-exp((log(y[i]) - mui)/new_sigma0[k1]))
    if (deni == 0)
    deni = 1e-06
    phi[i, k1] = new_pi0[k1] * deni
    sumi = sumi + phi[i, k1]
    }
    for (k1 in seq_len(Num.Comp)) {
    W[i, k1] = phi[i, k1] / sumi
    }
    }

    for (k1 in seq_len(Num.Comp)) {
    for (i in seq_len(Sample.Size)) {
    mui = sum(x[i,] * new_beta0[, k1] * acs[-1, k1]) +
    new_alpha0[k1] * acs[1, k1]
    predict[i, k1] = exp(mui)
    residual[i, k1] = y[i] - exp(mui)
    }
    }

    list(alpha = c(new_alpha0), beta = c(new_beta0),
    sigma = c(new_sigma0), pi = c(new_pi0), LogLik = loglike1,
    BIC = BIC, AIC = AIC, predict = c(predict),
    residual = c(residual), tau = c(W), MaxIter = MaxEMiter)
    }

    res <- fmrs.mle2(y = as.double(y),
    x = c(as.double(as.vector(unlist(x)))),
    delta = as.integer(delta),
    Lambda.Ridge = as.double(lambRidge),
    Num.Comp = as.integer(nComp),
    Num.Cov = as.integer(nCov),
    Sample.Size = as.integer(n),
    Num.iterationEM = as.integer(nIterEM),
    Num.iterationNR = 1, PortionNF = as.integer(NRpor),
    Initial.Intercept = as.double(c(coef0[, 1])),
    Initial.Coefficient = as.double(c(t(coef0[, -1]))),
    Initial.Dispersion = as.double(initDispersion),
    Initial.mixProp = as.double(initmixProp),
    conv.eps.em = as.double(convepsEM),
    actset = as.integer(activeset))

    print("The results are based on using survreg in survival package.")
    }

    } else {
    stop("The family of sub-distributions is not specified correctly.")
    }

    fit <- new("fmrsfit", y = y,
    delta = delta,
    x = x,
    nobs = n,
    ncov = nCov,
    ncomp = nComp,
    coefficients = array(rbind(res$alpha,
    matrix(res$beta, nrow = nCov, byrow = FALSE)), dim = c(nCov + 1, nComp),
    dimnames = list(xnames, comnames)),
    dispersion = array(res$sigma, dim = c(1, nComp),
    dimnames = list(NULL, comnames)),
    mixProp = array(res$pi, dim = c(1, nComp), dimnames = list(NULL, comnames)),
    logLik = res$LogLik,
    BIC = res$BIC,
    nIterEMconv = res$MaxIter,
    disFamily = disFamily,
    lambRidge = lambRidge,
    model = model,
    fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
    dim = c(n, nComp), dimnames = list(NULL,comnames)),
    residuals = array(matrix(res$residual, nrow = n, byrow = FALSE),
    dim = c(n, nComp), dimnames = list(NULL, comnames)),
    weights = array(matrix(res$tau, nrow = n, byrow = FALSE), dim = c(n, nComp),
    dimnames = list(NULL, comnames)),
    activeset = array(matrix(c(activeset), nrow = nCov + 1, byrow = FALSE),
    dim = c(nCov + 1, nComp), dimnames = list(xnames, comnames)),
    selectedset = array(matrix(1, nrow = nCov, ncol = nComp),
    dim = c(nCov, nComp), dimnames = list(xnames[-1], comnames))
    )
    return(fit)
    })

#' @rdname fmrs.tunsel-methods
#' @aliases fmrs.tunsel-method
setMethod(f = "fmrs.tunsel",
    definition = function(y, delta, x, nComp, disFamily = "lnorm",
    initCoeff, initDispersion, initmixProp,
    penFamily = "lasso", lambRidge = 0,
    nIterEM = 400, nIterNR = 2, conveps = 1e-08,
    convepsEM = 1e-08, convepsNR = 1e-08, NRpor = 2,
    gamMixPor = 1, activeset, lambMCP, lambSICA,
    cutpoint = 0.05, LambMin = 0.01, LambMax = 1.0,
    nLamb = 100) {
    if (missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
    if (missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
    if (missing(delta) & (disFamily != "norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

    if (missing(lambMCP))
    lambMCP = 2/(1 - max(cor(x)[cor(x) != 1]))

    if (missing(lambSICA))
    lambSICA = 5

    if ((nComp < 2)) {
    stop("An interger greater than 1 for the order of mixture model must be
    provided.")
    }

    if(nComp > 10){
    warning("A very high value is selected for nComp!")
    }
    nCov = dim(x)[2]
    n = length(y)
    if (missing(initCoeff))
    initCoeff = rnorm((nCov + 1) * nComp)
    if (missing(initDispersion))
    initDispersion = rep(1, nComp)
    if (missing(initmixProp))
    initmixProp = rep(1/nComp, nComp)
    if (missing(activeset))
    activeset = matrix(1, nrow = nCov + 1, ncol = nComp)
    if (any(activeset != 0 & activeset != 1))
    stop("activeset must be a matrix of dimention nCov+1 by nComp with
    only 0 and 1 values.")

    coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov + 1, byrow = TRUE)

    if (is.null(colnames(x))) {
    xnames <- c("Intercept", c(paste("X", seq_len(nCov), sep = ".")))
    } else {
    xnames <- c("Intercept", colnames(x))
    }
    comnames <- c(paste("Comp", seq_len(nComp), sep = "."))

    if (penFamily == "lasso")
    myPenaltyFamily = 1 else if (penFamily == "scad")
    myPenaltyFamily = 2 else if (penFamily == "mcp")
    myPenaltyFamily = 3 else if (penFamily == "sica")
    myPenaltyFamily = 4 else if (penFamily == "adplasso")
    myPenaltyFamily = 5 else if (penFamily == "hard")
    myPenaltyFamily = 6 else {
    stop("Penalty is not correctly specified.")
    }

    if (disFamily == "norm") {
    model = "FMR"
    delta = rep(1, n)
    res = .Call(FMR_Norm_CTun,
    as.double(y),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(c(coef0[,1])),
    as.double(c(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.integer(1),
    as.integer(myPenaltyFamily),
    as.double(gamMixPor),
    as.double(lambMCP),
    as.double(lambSICA),
    as.double(cutpoint),
    as.double(LambMin),
    as.double(LambMax),
    as.integer(nLamb)
    )
    } else if (disFamily == "lnorm") {
    model = "FMAFTR"
    logy = log(y)
    res = .Call(FMR_Norm_CTun,
    as.double(logy),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(c(coef0[,1])),
    as.double(c(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.integer(2),
    as.integer(myPenaltyFamily),
    as.double(gamMixPor),
    as.double(lambMCP),
    as.double(lambSICA),
    as.double(cutpoint),
    as.double(LambMin),
    as.double(LambMax),
    as.integer(nLamb)
    )
    } else if (disFamily == "weibull") {
    model = "FMAFTR"
    logy = log(y)
    res = .Call(FMR_Weibl_CTun,
    as.double(logy),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(c(coef0[,1])),
    as.double(c(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.integer(nIterNR),
    as.integer(myPenaltyFamily),
    as.double(gamMixPor),
    as.double(lambMCP),
    as.double(lambSICA),
    as.double(cutpoint),
    as.double(LambMin),
    as.double(LambMax),
    as.integer(nLamb),
    as.double(NRpor)
    )
    } else {
    stop("The family of sub-distributions is not specified correctly.")
    }

    lambdafit <- new("fmrstunpar", ncomp = nComp,
    ncov = nCov,
    lambPen = array(res$OptLam, dim = c(1, nComp),
    dimnames = c(list(NULL, c(paste("Comp", seq_len(nComp), sep = "."))))),
    lambRidge = lambRidge,
    disFamily = disFamily,
    penFamily = penFamily,
    MCPGam = lambMCP,
    SICAGam = lambSICA,
    model = model,
    activeset = array(matrix(c(activeset), nrow = nCov + 1, byrow = FALSE),
    dim = c(nCov + 1, nComp), dimnames = list(xnames, comnames))
    )
    return(lambdafit)
    })

#' @rdname fmrs.varsel-methods
#' @aliases fmrs.varsel-method
setMethod(f = "fmrs.varsel",
    definition = function(y, delta, x, nComp, disFamily = "lnorm",
    initCoeff, initDispersion, initmixProp,
    penFamily = "lasso", lambPen, lambRidge = 0,
    nIterEM = 2000, nIterNR = 2, conveps = 1e-08,
    convepsEM = 1e-08, convepsNR = 1e-08, NRpor = 2,
    gamMixPor = 1, activeset, lambMCP, lambSICA,
    cutpoint = 0.05) {
    if (missing(y) | !is.numeric(y))
    stop("A numeric response vector must be provided.")
    if (missing(x) | !is.numeric(x))
    stop("A numeric matrix for covariates must be provided.")
    if (missing(delta) & (disFamily != "norm"))
    stop("A censoring indicator vector with 0 or 1 values must be provided.")

    if (missing(lambMCP))
    lambMCP = 2/(1 - max(cor(x)[cor(x) != 1]))

    if (missing(lambSICA))
    lambSICA = 5

    if ((nComp < 2)) {
    stop("An interger greater than 1 for the order of mixture model must be
    provided.")
    }
    if(nComp > 10){
    warning("A very high value is selected for nComp!")
    }
    nCov = dim(x)[2]
    n = length(y)
    if (missing(initCoeff))
    initCoeff = rnorm((nCov + 1) * nComp)
    if (missing(initDispersion))
    initDispersion = rep(1, nComp)
    if (missing(initmixProp))
    initmixProp = rep(1/nComp, nComp)

    if (missing(activeset))
    activeset = matrix(1, nrow = nCov + 1, ncol = nComp)
    if (any(activeset != 0 & activeset != 1))
    stop("activeset must be a matrix of dimention nCov+1 by nComp with only 0
    and 1 values.")

    coef0 <- matrix(c(initCoeff), nrow = nComp, ncol = nCov+1, byrow = TRUE)

    if (is.null(colnames(x))) {
    xnames <- c("Intercept", c(paste("X", seq_len(nCov), sep = ".")))
    } else {
    xnames <- c("Intercept", colnames(x))
    }
    comnames <- c(paste("Comp", seq_len(nComp), sep = "."))

    if (penFamily == "lasso")
    myPenaltyFamily = 1 else if (penFamily == "scad")
    myPenaltyFamily = 2 else if (penFamily == "mcp")
    myPenaltyFamily = 3 else if (penFamily == "sica")
    myPenaltyFamily = 4 else if (penFamily == "adplasso")
    myPenaltyFamily = 5 else if (penFamily == "hard")
    myPenaltyFamily = 6 else {
    stop("Penalty is not correctly specified.")
    }

    if (disFamily == "norm")
    {
    model = "FMR"
    delta = rep(1, n)

    res = .Call(FMR_Norm_MPLE,
    as.double(y),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(c(coef0[, 1])),
    as.double(c(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.double(convepsEM),
    as.integer(1),
    as.integer(myPenaltyFamily),
    as.double(lambPen),
    as.double(gamMixPor),
    as.double(lambMCP),
    as.double(lambSICA),
    as.double(cutpoint)
    )
    } else if (disFamily == "lnorm")
    {
    model = "FMAFTR"
    logy = log(y)
    res = .Call(FMR_Norm_MPLE,
    as.double(logy),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(c(coef0[, 1])),
    as.double(c(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.double(convepsEM),
    as.integer(2),
    as.integer(myPenaltyFamily),
    as.double(lambPen),
    as.double(gamMixPor),
    as.double(lambMCP),
    as.double(lambSICA),
    as.double(cutpoint)
    )
    } else if (disFamily == "weibull")
    {
    model = "FMAFTR"
    logy = log(y)
    res = .Call(FMR_Weibl_MPLE,
    as.double(logy),
    as.double(as.vector(unlist(x))),
    as.integer(nComp),
    as.integer(nCov),
    as.integer(n),
    as.integer(delta),
    as.double(c(coef0[, 1])),
    as.double(c(t(coef0[, -1]))),
    as.double(initDispersion),
    as.double(initmixProp),
    as.integer(activeset),
    as.double(lambRidge),
    as.integer(nIterEM),
    as.double(conveps),
    as.double(convepsEM),
    as.integer(2),
    as.integer(myPenaltyFamily),
    as.double(lambPen),
    as.double(gamMixPor),
    as.double(lambMCP),
    as.double(lambSICA),
    as.double(cutpoint),
    as.double(NRpor)
    )
    } else {
    stop("The family of sub-distributions is not specified correctly.")
    }

    fit <- new("fmrsfit", y = y, delta = delta, x = x, nobs = n, ncov = nCov,
    ncomp = nComp,
    coefficients = array(rbind(res$alpha,
    matrix(res$beta,nrow = nCov, byrow = FALSE)), dim = c(nCov + 1, nComp),
    dimnames = list(xnames, comnames)),
    dispersion = array(res$sigma, dim = c(1, nComp),
    dimnames = list(NULL, comnames)),
    mixProp = array(res$pi, dim = c(1, nComp), dimnames = list(NULL, comnames)),
    logLik = res$LogLik,
    BIC = res$BIC,
    nIterEMconv = res$MaxIter,
    disFamily = disFamily,
    penFamily = penFamily,
    lambPen = array(lambPen, dim = c(1, nComp), dimnames = list(NULL,comnames)),
    MCPGam = lambMCP,
    SICAGam = lambSICA,
    model = model,
    fitted = array(matrix(res$predict, nrow = n, byrow = FALSE),
    dim = c(n, nComp), dimnames = list(NULL, comnames)),
    residuals = array(matrix(res$residual,nrow = n, byrow = FALSE),
    dim = c(n, nComp), dimnames = list(NULL, comnames)),
    weights = array(matrix(res$tau, nrow = n, byrow = FALSE),
    dim = c(n, nComp), dimnames = list(NULL, comnames)),
    activeset = array(matrix(c(activeset), nrow = nCov + 1, byrow = FALSE),
    dim = c(nCov + 1, nComp), dimnames = list(xnames,comnames)),
    selectedset = array(matrix(c(res$selection), nrow = nCov, byrow = FALSE),
    dim = c(nCov, nComp), dimnames = list(xnames[-1],comnames))
    )
    return(fit)
    })

#' @rdname fmrs.gendata-methods
#' @aliases fmrs.gendata-method
setMethod(f = "fmrs.gendata",
    definition = function(nObs, nComp, nCov, coeff, dispersion, mixProp,
    rho, umax, disFamily = "lnorm") {
    if (missing(disFamily))
    disFamily = "lnorm"

    if (sum(mixProp) != 1)
    stop("The sum of mixing proportions must be 1.")

    if (sum(dispersion <= 0) != 0)
    stop("Dispersion parameters cannot be zero or negative.")

    if (rho > 1 | rho < -1)
    stop("The correlation cannot be less than -1 or greater thatn 1.")

    if ((nComp*(nCov+1)) != length(coeff))
    stop("The length of coeff is different from nComp *(nCov+1).")

    if (nComp != length(dispersion))
    stop("The length of dispersian is different from nComp.")

    if (nComp != length(mixProp))
    stop("The length of mixProp is different from nComp.")

    if ((disFamily != "norm") & (umax<=0))
    stop("Provide a positive value for umax.")

    if(nComp > 10){
    warning("A very high value is selected for nComp!")
    }

    mu <- rep(0, nCov)
    Sigma <- diag(nCov)

    for (i in seq_len(nCov)) {
    for (j in seq_len(nCov)) {
    Sigma[i, j] <- rho^abs(i - j)
    }
    }

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
    cX = drop(X) else cX = t(X)
    cX <- scale(cX)
    colnames(cX) <- paste("X", seq_len(nCov), sep = ".")

    coef0 <- matrix(coeff, nrow = nComp, ncol = nCov + 1, byrow = TRUE)
    mixProp0 <- cumsum(mixProp)

    yobs <- c()
    c <- rep()
    dlt <- c()
    u <- c()
    tobs <- c()

    if (disFamily == "lnorm") {
    for (i in seq_len(nObs)) {
    eps <- rnorm(1)
    u1 <- runif(1)
    k = length(which(mixProp0 <= u1)) + 1
    u[i] = k
    yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * eps

    c[i] <- log(runif(1, 0, umax))
    tobs[i] <- exp(min(yobs[i], c[i]))
    dlt[i] <- (yobs[i] < c[i]) * 1
    }
    } else if (disFamily == "norm") {
    for (i in seq_len(nObs)) {
    eps <- rnorm(1)
    u1 <- runif(1)
    k = length(which(mixProp0 <= u1)) + 1
    u[i] = k
    yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * eps
    tobs[i] <- yobs[i]
    dlt[i] <- 1
    }
    } else if (disFamily == "weibull") {
    for (i in seq_len(nObs)) {
    ext <- log(rexp(1))
    u1 <- runif(1)
    k = length(which(mixProp0 <= u1)) + 1
    yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * ext

    c[i] <- log(runif(1, 0, umax))
    tobs[i] <- exp(min(yobs[i], c[i]))
    dlt[i] <- (yobs[i] < c[i]) * 1
    }
    } else {
    stop("The family of sub-distributions are not specified correctly.")
    }
    return(list(y = tobs, delta = dlt, x = cX, disFamily = disFamily))
    })
