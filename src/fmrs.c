#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <stdio.h>

/* ************************************************************************ */
/* .Call calls */
extern SEXP FMR_Norm_CTun(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FMR_Norm_MLE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FMR_Norm_MPLE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FMR_Weibl_CTun(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FMR_Weibl_MLE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FMR_Weibl_MPLE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"FMR_Norm_CTun",  (DL_FUNC) &FMR_Norm_CTun,  23},
    {"FMR_Norm_MLE",   (DL_FUNC) &FMR_Norm_MLE,   16},
    {"FMR_Norm_MPLE",  (DL_FUNC) &FMR_Norm_MPLE,  22},
    {"FMR_Weibl_CTun", (DL_FUNC) &FMR_Weibl_CTun, 24},
    {"FMR_Weibl_MLE",  (DL_FUNC) &FMR_Weibl_MLE,  17},
    {"FMR_Weibl_MPLE", (DL_FUNC) &FMR_Weibl_MPLE, 23},
    {NULL, NULL, 0}
};

void R_init_fmrs(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

/* ************************************************************************ */
/* ************************************************************************ */
void backsub(int n, double *a, double *y)
{
    int i, k, j;
    double atemp[n + 1][n + 1];

    for (i=0; i < (n + 1); i++) {
    for (j=0; j < (n + 1); j++) {
    atemp[i][j] = a[j + i * (n + 1)];
    }
    }

    y[n-1] = atemp[n - 1][n] / atemp[n - 1][n - 1];

    for(i = n - 2; i >= 0; i--) {
    y[i] = atemp[i][n];
    for(k = n - 1; k > i; k--)
    y[i] = y[i] - atemp[i][k] * y[k];
    y[i] = y[i] / atemp[i][i];
    }
}

/* ************************************************************************ */
void maxabs(int row, int n, double *a, int *checksin, int *maxrow)
{
    double max=0;
    int loc_maxrow, loc_checksin;
    int i, j;
    loc_maxrow = row;
    loc_checksin = 1;
    double atemp[n + 1][n + 1];

    for (i=0; i < (n + 1); i++) {
    for (j=0; j < (n + 1); j++) {
    atemp[i][j] = a[j + i * (n + 1)];
    }
    }

    for (i = row; i < n; i++) {
    if ( fabs(atemp[i][row]) > max) {
    max = fabs(atemp[i][row]);
    loc_maxrow = i;
    }

    if (max == 0) {
    loc_checksin=0;
    break;
    }
    }
    maxrow[0]   = loc_maxrow;
    checksin[0] = loc_checksin;
}

/* ************************************************************************ */
void gauss(int n, double *a, int *checksin)
{
    int i, j, k;
    int k_vec[1];
    double temp;
    double atemp[n + 1][n + 1];
    k=0;

    for (i=0; i < (n + 1); i++) {
    for (j=0; j < (n + 1); j++) {
    atemp[i][j] = a[j + i * (n + 1)];
    }
    }

    for (i=0; i < n; i++) {
    maxabs(i, n, a, checksin, k_vec);
    k = k_vec[0];

    if (checksin[0] == 0) {
    break;
    }

    if (k != i)
    for (j = i; j < n + 1; j++) {
    temp = atemp[i][j];
    atemp[i][j] = atemp[k][j];
    atemp[k][j] = temp;
    }

    for (j = i + 1; j < n; j++) {
    temp = atemp[j][i] / atemp[i][i];
    for (k = i; k < n + 1; k++)
    atemp[j][k] = atemp[j][k] - 1.0 * temp * atemp[i][k];
    }
    }

    for (i=0; i < (n + 1); i++) {
    for (j=0; j < (n + 1); j++) {
    a[j + i * (n + 1)] = atemp[i][j];
    }
    }
}

/* ************************************************************************ */
void sol(int Nelem, double *IS, double *solution, int *checksin )
{
    gauss(Nelem, IS, checksin);
    backsub(Nelem, IS, solution);
}

/* ************************************************************************ */
double minimum(double *vector1, int D1)
{
    int i;
    double min1;
    min1 = fabs(vector1[0]);

    for(i = 1; i < D1; i++)
    if(fabs(vector1[i]) <= min1)
    min1 = fabs(vector1[i]);

    return min1;
}
/* ************************************************************************ */
double maximum(double *vector1, int D1)
{
    int i;
    double max1;
    max1 = fabs(vector1[0]);

    for(i=0; i < D1; i++)
    if(fabs(vector1[i]) >= max1)
    max1 = fabs(vector1[i]);
    return max1;
}

/* ************************************************************************ */
void sicapen(double lam1, double *moshtagh, double *regcoef, int D1,
     double SicaTun)
{
    int j;

    for(j=0; j < D1; j++)
    moshtagh[j] = lam1 * SicaTun * (SicaTun + 1) / pow(SicaTun +
    fabs(regcoef[j]), 2);
}

/* ************************************************************************ */
void mcppen(double lam1, double *moshtagh, double *regcoef, int D1,
    double McpTun)
{
    int j;

    for(j=0; j < D1; j++){
    if(fabs(regcoef[j]) <= (lam1 * McpTun))
    moshtagh[j] = lam1 * (1 - fabs(regcoef[j]) / (lam1 * McpTun));
    else
    moshtagh[j]= 0.0;
    }

}

/* ************************************************************************ */
void scadpen(double lam1, double *moshtagh, double *regcoef, int D1)
{
    int j;
    double a1;
    a1=3.7;

    for(j=0; j < D1; j++){
    if(fabs(regcoef[j]) <= lam1)
    moshtagh[j] = lam1;
    else if((a1 * lam1 - fabs(regcoef[j])) <= 0)
    moshtagh[j]=0;
    else
    moshtagh[j]= (a1 * lam1 - fabs(regcoef[j])) / (a1 - 1);
    }
}

/* ************************************************************************ */
void hardpen(double lam1, double *moshtagh, double *regcoef, int D1)
{
    int j;

    for(j=0; j < D1; j++){
    if(fabs(regcoef[j]) <= lam1)
    moshtagh[j]= (-2) * (fabs(regcoef[j]) - lam1);
    else
    moshtagh[j]=0;
    }
}

const char* names[] = { "alpha", "beta", "sigma", "pi", "LogLik", "BIC",
    "AIC", "MaxIter", "tau", "predict", "residual", "" };

int accessAcsArr(SEXP acs, const int first_index, const int second_index,
    const int D)
{
    return INTEGER(acs)[(D + 1) * second_index + first_index];
}

/* *************************** MLE Normal and Log-Normal ******************* */
SEXP FMR_Norm_MLE(SEXP myY,
      SEXP myX,
      SEXP myK,
      SEXP myD,
      SEXP myN,
      SEXP mydelta,
      SEXP myAlpha,
      SEXP myBeta,
      SEXP mySigma,
      SEXP myPi,
      SEXP myacs,
      SEXP myridge,
      SEXP myEMiter,
      SEXP myeps,
      SEXP myepsC,
      SEXP myNorm
)
{
    int i;
    int j;
    int k1;
    int l;
    int niter1;
    int check1[1];

    int N = asInteger(myN);
    int D = asInteger(myD);
    int EM_Miter = asInteger(myEMiter);
    int K = asInteger(myK);

    int OD[K];
    for (k1 = 0; k1 < K; k1++) {
    OD[k1] = 0;
    for (j = 0; j < D + 1; j++)
    OD[k1] += accessAcsArr(myacs, j, k1, D);
    }

    double eps_conv = asReal(myepsC);
    double eps = asReal(myeps);

    double ridge = asReal(myridge);
    double one_X[N][D + 1][K];
    double sumi;
    double mui;
    double deni;

    double W[N][K];
    double phi[N][K];

    double beta_ini[D][K];
    double beta_old[D][K];
    double beta_new[D][K];

    double alpha_ini[K];
    double alpha_old[K];
    double alpha_new[K];

    double sigma_ini[K];
    double sigma_old[K];
    double sigma_new[K];

    double pi_ini[K];
    double pi_old[K];
    double pi_new[K];

    double w_s[N][K];
    double Aw[N][K];
    double V[N][K];

    double oneXTWY[D + 1];
    double oneXTWX[D + 1][D + 1];

    double sumwi[K];
    double sumi5[K];
    double sumi3[K];

    double SumDif;

    double oneComMat[D + 2][D + 2];
    double OneSolv[D + 2];
    double oneComMatVec[(D + 2) * (D + 2)];
    double loglike1;

    double sumi1, sumi2;

    for (k1 = 0; k1 < K; k1++) {
    for (j = 0; j < D; j++) {
    beta_new[j][k1] = beta_old[j][k1] =
    beta_ini[j][k1] = REAL(myBeta)[D * k1 + j] *
    accessAcsArr(myacs, j + 1, k1, D);
    }
    alpha_new[k1] = alpha_old[k1] = alpha_ini[k1] =
    REAL(myAlpha)[k1] * accessAcsArr(myacs, 0, k1, D);
    sigma_new[k1] = sigma_old[k1] = sigma_ini[k1] =
    REAL(mySigma)[k1];
    pi_new[k1] = pi_old[k1] = pi_ini[k1] = REAL(myPi)[k1];
    }

    for(k1=0; k1<K; k1++){
    l = 0;
    if (accessAcsArr(myacs, 0, k1, D) == 1) {
    for (i = 0; i < N; i++) {
    one_X[i][l][k1] = 1.0;
    }
    l++;
    }

    for (j = 0; j < D; j++) {
    if (accessAcsArr(myacs, j + 1, k1, D) == 1) {
    for (i = 0; i < N; i++) {
    one_X[i][l][k1] = REAL(myX)[N * j + i];
    }
    l++;
    }
    }
    }

    niter1 = 0;
    while (niter1 < EM_Miter) {
    /****Beginning of each iteration****/
    /*******The E-step of the EM********/

    for (k1 = 0; k1 < K; k1++) {
    sumwi[k1] = 0.0;
    sumi5[k1] = 0.0;
    sumi3[k1] = 0.0;
    }

    for (i = 0; i < N; i++) {
    sumi = 0.0;
    for (k1 = 0; k1 < K; k1++) {
    mui = 0.0;
    for (j = 0; j < D; j++)
    mui += REAL(myX)[N * j + i] * beta_old[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    deni = pow(eps + dnorm(REAL(myY) [i] - mui, 0, sigma_old[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_old[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    phi[i][k1] = pi_old[k1] * deni;
    sumi += phi[i][k1];
    }
    for (k1 = 0; k1 < K; k1++) {
    W[i][k1] = phi[i][k1] / sumi;
    sumwi[k1] += W[i][k1];
    }
    }

    /**********End of the E-step*******/

    for (i = 0; i < N; i++) {
    for (k1 = 0; k1 < K; k1++) {
    mui = 0.0;
    for (j = 0; j < D; j++)
    mui += REAL(myX)[N * j + i] * beta_old[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    w_s[i][k1] = (REAL(myY)[i] - mui) / sigma_old[k1];
    Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 -
    pnorm(w_s[i][k1], 0, 1, 1, 0)));
    V[i][k1] = INTEGER(mydelta)[i] * REAL(myY)[i] +
    (1 - INTEGER(mydelta)[i]) * (mui + sigma_old[k1] * Aw[i][k1]);
    }
    }

    /*********M-step of the EM*********/
    for (k1 = 0; k1 < K; k1++) {

    /******Constructing the weighted vector XTWY***/

    for (j = 0; j < OD[k1]; j++) {
    sumi1 = 0.0;
    for (i = 0; i < N; i++)
    sumi1 += one_X[i][j][k1] * W[i][k1] * V[i][k1];
    oneXTWY[j] = sumi1;
    }

    /*****Constructing the weighted matrix XTWX***/

    for (i = 0; i < (OD[k1]); i++) {
    for (j = i; j < (OD[k1]); j++) {
    sumi2 = 0.0;
    for (l = 0; l < N; l++)
    sumi2 += one_X[l][i][k1] * W[l][k1] * one_X[l][j][k1];
    if (i == j)
    oneXTWX[i][j] = sumi2 + ridge * log(N);
    else
    oneXTWX[j][i] = oneXTWX[i][j] = sumi2;
    }
    }

    if (accessAcsArr(myacs, 0, k1, D) == 1)
    oneXTWX[0][0] = oneXTWX[0][0] - ridge * log(N);

    /***In a system Ax=b, adding b to A as its last column**/

    for (i = 0; i < OD[k1]; i++) {
    for (j = 0; j < (OD[k1] + 1); j++) {
    if (j != OD[k1])
    oneComMat[i][j] = oneXTWX[i][j];
    else
    oneComMat[i][j] = oneXTWY[i];
    }
    }

    for (i = 0; i < (OD[k1] + 1); i++) {
    for (j = 0; j < (OD[k1] + 1); j++) {
    oneComMatVec[j + i * (OD[k1] + 1)] = oneComMat[i][j];
    }
    }

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/
    /**************************************************************/

    sol(OD[k1], oneComMatVec, OneSolv, check1);

    l = 0;
    if (accessAcsArr(myacs, l, k1, D) == 1)
    alpha_new[k1] = OneSolv[l++];
    for (j = 0; j < D; j++) {
    if (accessAcsArr(myacs, j + 1, k1, D) == 1) {
    beta_new[j][k1] = OneSolv[l++];
    }
    }
    }

    for (i = 0; i < N; i++) {
    for (k1 = 0; k1 < K; k1++) {
    mui = 0.0;
    for (j = 0; j < D; j++)
    mui += REAL(myX)[N * j + i] * beta_new[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    sumi5[k1] += W[i][k1] * pow(V[i][k1] - mui, 2);
    sumi3[k1] += W[i][k1] * (INTEGER(mydelta)[i] + (1 - INTEGER(mydelta)[i]) *
    (Aw[i][k1] * (Aw[i][k1] - w_s[i][k1])));
    }
    }

    for (k1 = 0; k1 < K; k1++) {
    sigma_new[k1] = sqrt((sumi5[k1]) / (sumi3[k1]));
    }

    for (k1 = 0; k1 < K; k1++) {
    pi_new[k1] = sumwi[k1] / N;
    }

    /*****End of the M-step of the EM*****/
    SumDif = 0.0;
    niter1++;
    for (k1 = 0; k1 < K; k1++) {
    for (j = 0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_old[j][k1], 2) *
    accessAcsArr(myacs, j + 1, k1, D);
    SumDif += pow(alpha_new[k1] - alpha_old[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);
    SumDif += pow(pi_new[k1] - pi_old[k1], 2);
    SumDif += pow(sigma_new[k1] - sigma_old[k1], 2);
    }

    if (SumDif <= eps_conv)
    break;

    for (k1 = 0; k1 < K; k1++) {
    for (j = 0; j < D; j++)
    beta_old[j][k1] = beta_new[j][k1];
    alpha_old[k1] = alpha_new[k1];
    pi_old[k1] = pi_new[k1];
    sigma_old[k1] = sigma_new[k1];
    }

    /*******End of each iteration*******/
    }

    loglike1 = 0.0;

    for (i = 0; i < N; i++) {
    sumi = 0.0;
    for (k1 = 0; k1 < K; k1++) {
    mui = 0.0;
    for (j = 0; j < D; j++)
    mui += REAL(myX)[N * j + i] * beta_new[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    deni = pow(eps + dnorm(REAL(myY)[i] - mui, 0, sigma_new[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_new[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    phi[i][k1] = pi_new[k1] * deni;
    sumi += phi[i][k1];
    }
    loglike1 += log(sumi);
    for (k1 = 0; k1 < K; k1++) {
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    int SumPar = 0;
    for (k1 = 0; k1 < K; k1++) {
    SumPar += OD[k1];
    }
    SumPar += 2 * K - 1;

    double BIC = 0.0;
    double AIC = 0.0;
    BIC = -(2 * loglike1) + SumPar * log(N);
    AIC = -(2 * loglike1) + 2 * SumPar;

    SEXP alpha = PROTECT(allocVector(REALSXP, K));
    SEXP beta = PROTECT(allocVector(REALSXP, K * D));
    SEXP sigma = PROTECT(allocVector(REALSXP, K));
    SEXP pi = PROTECT(allocVector(REALSXP, K));
    SEXP predict = PROTECT(allocVector(REALSXP, K * N));
    SEXP residual = PROTECT(allocVector(REALSXP, K * N));
    SEXP tau = PROTECT(allocVector(REALSXP, K * N));


    for (k1 = 0; k1 < K; k1++) {
    for (j = 0; j < D; j++)
    REAL(beta)[k1 * D + j] = beta_new[j][k1];
    REAL(alpha)[k1] = alpha_new[k1];
    REAL(sigma)[k1] = sigma_new[k1];
    REAL(pi)[k1] = pi_new[k1];
    }

    /**********End of the E-step*******/
    for (k1 = 0; k1 < K; k1++) {
    for (i = 0; i < N; i++) {
    mui = 0.0;
    for (j = 0; j < D; j++)
    mui += REAL(myX)[N * j + i] * beta_new[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);

    if (asInteger(myNorm) == 1) {
    REAL(predict)[k1 * N + i] = mui;
    REAL(residual)[k1 * N + i] = REAL(myY)[i] - mui;
    }
    else {
    REAL(predict)[k1 * N + i] = exp(mui);
    REAL(residual)[k1 * N + i] = exp(REAL(myY)[i]) - exp(mui);
    }
    REAL(tau)[k1 * N + i] = W[i][k1];
    }
    }

    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, alpha);
    SET_VECTOR_ELT(res, 1, beta);
    SET_VECTOR_ELT(res, 2, sigma);
    SET_VECTOR_ELT(res, 3, pi);
    SET_VECTOR_ELT(res, 4, ScalarReal(loglike1));
    SET_VECTOR_ELT(res, 5, ScalarReal(BIC));
    SET_VECTOR_ELT(res, 6, ScalarReal(AIC));
    SET_VECTOR_ELT(res, 7, ScalarInteger(niter1));
    SET_VECTOR_ELT(res, 8, tau);
    SET_VECTOR_ELT(res, 9, predict);
    SET_VECTOR_ELT(res, 10, residual);
    UNPROTECT(8);
    return res;
}

/* ********************* Variable Selection Normal and Log-Normal *********** */
SEXP FMR_Norm_MPLE(SEXP myY,
       SEXP myX,
       SEXP myK,
       SEXP myD,
       SEXP myN,
       SEXP mydelta,
       SEXP myAlpha,
       SEXP myBeta,
       SEXP mySigma,
       SEXP myPi,
       SEXP myacs,
       SEXP myridge,
       SEXP myEMiter,
       SEXP myeps,
       SEXP myepsC,
       SEXP myNorm,
       SEXP myPenalty,
       SEXP mylambda,
       SEXP myPiPor,
       SEXP myMcpTun,
       SEXP mySicaTun,
       SEXP myTol
)
{
    int i;
    int j;
    int k1;
    int l;
    int niter1;
    int check1[1];
    int SUM1;

    int N = asInteger(myN);
    int D = asInteger(myD);
    int EM_Miter = asReal(myEMiter);
    int K = asInteger(myK);
    int jar = asInteger(myPenalty);
    double McpTun = asReal(myMcpTun);
    double SicaTun = asReal(mySicaTun);

    double Tol = asReal(myTol);

    int OD[K];
    for(k1=0; k1<K; k1++){
    OD[k1]=0;
    for(j=0; j<D+1; j++){
    OD[k1] += accessAcsArr(myacs, j, k1, D);
    }
    }

    double ridge = asReal(myridge);
    double eps = asReal(myeps);
    double eps_conv = asReal(myepsC);
    double GamMP = asReal(myPiPor);

    double one_X[N][D+1][K];

    double optlam[K];

    double beta_ini[D][K];
    double beta_old[D][K];
    double beta_new[D][K];

    double alpha_ini[K];
    double alpha_old[K];
    double alpha_new[K];

    double sigma_ini[K];
    double sigma_old[K];
    double sigma_new[K];

    double pi_ini[K];
    double pi_old[K];
    double pi_new[K];

    double W[N][K];
    double phi[N][K];

    double w_s[N][K];
    double Aw[N][K];
    double V[N][K];

    double eps1[K];

    double vecder[D];
    double vecSig[D+1];
    double En[K];

    double oneXTWY[D+1];
    double oneXTWX[D+1][D+1];
    double selection[D][K];

    double oneComMat[D+2][D+2];
    double OneSolv[D+2];
    double oneComMatVec[(D+2) * (D+2)];

    double sumwi[K];
    double sumi4[K];
    double sumi5[K];
    double sumi;
    double mui;
    double deni;

    char convg1 = 'n';

    double SumDif;
    double loglike1;
    double holdveccov[D];

    for (k1=0; k1 < K; k1++) {
    optlam[k1]= REAL(mylambda)[k1];
    }

    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++){
    beta_new[j][k1]= beta_old[j][k1] = beta_ini[j][k1] =
    REAL(myBeta)[D*k1+j] *accessAcsArr(myacs, j + 1, k1, D) ;
    }
    alpha_new[k1] = alpha_old[k1] = alpha_ini[k1] =
    REAL(myAlpha)[k1] * accessAcsArr(myacs, 0, k1, D);
    sigma_new[k1] = sigma_old[k1] = sigma_ini[k1] =
    REAL(mySigma)[k1];
    pi_new[k1] = pi_old[k1] = pi_ini[k1] =
    REAL(myPi)[k1];
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    if(accessAcsArr(myacs, 0, k1, D)  == 1){
    for (i=0; i < N; i++){
    one_X[i][l][k1] = 1.0;
    }
    l++;
    }
    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j + 1, k1, D)  == 1){
    for (i=0; i < N; i++){
    one_X[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    for (j=0; j < D; j++) {
    holdveccov[j] = beta_ini[j][k1];
    }
    if((jar == 1) || (jar == 2))
    eps1[k1] = eps * minimum(holdveccov, D) / (2 * N * optlam[k1]);
    else
    eps1[k1] = eps * minimum(holdveccov, D) / (4 * N * optlam[k1]);
    }

    niter1=0;
    convg1 = 'n';

    while((convg1 != 'y') && (niter1 < EM_Miter)){
    /******Beginning of each iteration******/
    /****E-step of the EM algorithm************/
    for(k1=0; k1 < K;  k1++){
    sumwi[k1] = 0.0;
    }

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K;  k1++){
    mui=0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D) ;
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    deni = pow(eps + dnorm(REAL(myY)[i] - mui, 0, sigma_old[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_old[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    phi[i][k1] = pi_old[k1] * deni;
    sumi += phi[i][k1];
    }
    for(k1=0; k1 < K;  k1++){
    W[i][k1] = phi[i][k1] / sumi;
    sumwi[k1] += W[i][k1];
    }
    }

    for (k1=0; k1 < K; k1++) {
    pi_new[k1] = sumwi[k1] / N;
    }

    /*****End of the E-step of the EM algorithm*****/

    for(i=0; i < N; i++){
    for(k1=0; k1 < K;  k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D) ;
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    w_s[i][k1] = (REAL(myY)[i] - mui) / sigma_old[k1];
    Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 -
    pnorm(w_s[i][k1], 0, 1, 1, 0)));
    V[i][k1] = INTEGER(mydelta)[i] * REAL(myY)[i] +
    (1- INTEGER(mydelta)[i]) * (mui + sigma_old[k1] * Aw[i][k1]);
    }
    }

    /****M-step of the EM algorithm*************/
    for(k1=0; k1 < K;  k1++){

    En[k1] = N * pow(pi_old[k1], GamMP);

    for (j=0; j < D; j++) {
    holdveccov[j] = beta_old[j][k1];
    }

    if(jar == 1)
    for(j=0; j < D; j++)
    vecder[j] = optlam[k1];
    else if (jar == 2)
    scadpen(optlam[k1], vecder, holdveccov, D);
    else if (jar == 3)
    mcppen(optlam[k1], vecder, holdveccov, D, McpTun);
    else if (jar == 4)
    sicapen(optlam[k1], vecder, holdveccov, D, SicaTun);
    else if (jar == 5){
    for(j=0; j < D; j++)
    vecder[j] = optlam[k1] / (fabs(beta_ini[j][k1]) + eps);
    }
    else
    hardpen(optlam[k1], vecder, holdveccov, D);


    l=0;
    if(accessAcsArr(myacs, 0, k1, D)  == 1)
    vecSig[l++] = 0.0;
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j + 1, k1, D)  == 1)
    vecSig[l++] = vecder[j] / (fabs(beta_old[j][k1]) + eps1[k1]);
    }

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < OD[k1]; j++){
    sumi = 0.0;
    for(l=0; l < N; l++)
    sumi +=one_X[l][i][k1] * W[l][k1] * one_X[l][j][k1];
    if(i != j)
    oneXTWX[i][j] = sumi;
    else
    oneXTWX[i][j] = sumi + En[k1] * vecSig[j] +
    ridge * log(N);
    }
    }

    if(accessAcsArr(myacs, 0, k1, D) == 1)
    oneXTWX[0][0] = oneXTWX[0][0] - (En[k1] * vecSig[0]) -
    ridge * log(N);

    /******Constructing the weigthed vector XTWY***/

    for(j=0; j < OD[k1]; j++){
    sumi = 0.0;
    for(i=0; i < N; i++)
    sumi += one_X[i][j][k1] * W[i][k1] * V[i][k1];
    oneXTWY[j] = sumi;
    }

    /***In a system Ax=b, adding b to A as its last column**/

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < (OD[k1] + 1); j++){
    if(j != OD[k1])
    oneComMat[i][j] = oneXTWX[i][j];
    else
    oneComMat[i][j] = oneXTWY[i];
    }
    }

    for (i=0; i < (OD[k1] + 1); i++) {
    for (j=0; j < (OD[k1] + 1); j++) {
    oneComMatVec[j+i*(OD[k1]+1)] = oneComMat[i][j];
    }
    }

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/
    /**************************************************************/
    sol(OD[k1], oneComMatVec, OneSolv, check1);

    l=0;
    if(accessAcsArr(myacs, 0, k1, D)  == 1)
    alpha_new[k1] = OneSolv[l++];
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j + 1, k1, D)  == 1){
    beta_new[j][k1]  = OneSolv[l++];
    }
    }

    }

    for (k1=0; k1 < K; k1++){
    sumi5[k1] = 0.0;
    sumi4[k1] = 0.0;
    }

    for(i=0; i < N; i++){
    for(k1=0; k1 < K;  k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D) ;
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    sumi5[k1] += W[i][k1] * pow(V[i][k1] - mui, 2);
    }
    }
    for (k1=0; k1 < K; k1++) {
    for (i=0; i < N; i++) {
    sumi4[k1] += W[i][k1] * (INTEGER(mydelta)[i] +
    (1 - INTEGER(mydelta)[i]) * (Aw[i][k1] * (Aw[i][k1] - w_s[i][k1])));
    }
    }

    for (k1=0; k1 < K; k1++) {
    sigma_new[k1] = sqrt((sumi5[k1] / sumi4[k1]));
    }

    SumDif = 0.0;
    for(k1=0; k1 < K;  k1++){
    for(j=0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_old[j][k1], 2) *
    accessAcsArr(myacs, j + 1, k1, D);
    SumDif += pow(alpha_new[k1] - alpha_old[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);
    SumDif += pow(pi_new[k1] - pi_old[k1], 2);
    SumDif += pow(sigma_new[k1] - sigma_old[k1], 2);
    }

    if (SumDif <= eps_conv)
    convg1='y';

    niter1++;

    for(k1=0; k1 < K;  k1++){
    alpha_old[k1] = alpha_new[k1];
    for(j=0; j < D; j++)
    beta_old[j][k1] = beta_new[j][k1];
    pi_old[k1] = pi_new[k1];
    sigma_old[k1] = sigma_new[k1];
    }
    }

    //*********************************************************************
    //*Storing the estimates of the regression coefficients
    //*in a global variable called "Betahat".
    //*Selecting the finial model and storing
    //*********************************************************************

    for(k1=0; k1 < K;  k1++){
    for(j=0; j < D; j++)
    selection[j][k1]=0;
    }

    for(k1=0; k1 < K;  k1++){
    for(j=0; j < D; j++){
    if(fabs(beta_new[j][k1]) <= Tol){
    selection[j][k1]=0;
    beta_new[j][k1] = 0.0;
    }
    else
    selection[j][k1] = 1;
    }
    }


    for(i=0; i < N; i++){
    for(k1=0; k1 < K;  k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    selection[j][k1] * accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    w_s[i][k1] = (REAL(myY)[i] - mui) / sigma_new[k1];
    Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 -
    pnorm(w_s[i][k1], 0, 1, 1, 0)));
    V[i][k1] = INTEGER(mydelta)[i] * REAL(myY)[i] +
    (1 - INTEGER(mydelta)[i]) * (mui + sigma_new[k1] * Aw[i][k1]);
    }
    }

    loglike1 = 0.0;

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K;  k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    selection[j][k1] * accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    deni = pow(eps +  dnorm(REAL(myY)[i] - mui, 0, sigma_new[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_new[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    phi[i][k1] = pi_new[k1] * deni;
    sumi += phi[i][k1];
    }
    loglike1 += log(sumi);
    }

    SUM1=0;
    for(k1=0; k1 < K;  k1++){
    SUM1 += accessAcsArr(myacs, 0, k1, D);
    for(j=0; j < D; j++){
    SUM1 += selection[j][k1];
    }
    }
    SUM1 +=  2*K-1;
    double BIC = 0.0;
    double AIC = 0.0;

    BIC = -(2 * loglike1) + SUM1 * log(N);
    AIC = -(2 * loglike1) + 2 * SUM1;

    SEXP alpha = PROTECT(allocVector(REALSXP, K));
    SEXP beta = PROTECT(allocVector(REALSXP, K*D));
    SEXP sigma = PROTECT(allocVector(REALSXP, K));
    SEXP pi = PROTECT(allocVector(REALSXP, K));
    SEXP predict = PROTECT(allocVector(REALSXP, K*N));
    SEXP residual = PROTECT(allocVector(REALSXP, K*N));
    SEXP tau = PROTECT(allocVector(REALSXP, K*N));
    SEXP SeL = PROTECT(allocVector(INTSXP, K*D));

    for (k1=0; k1 < K; k1++){
    for (j=0;  j < D;  j++)
    INTEGER(SeL)[k1 * D + j] = selection[j][k1];
    }

    for (k1=0; k1 < K; k1++){
    for (j=0;  j < D;  j++)
    REAL(beta)[k1 * D + j] = beta_new[j][k1];
    REAL(alpha)[k1] = alpha_new[k1];
    REAL(sigma)[k1] = sigma_new[k1];
    REAL(pi)[k1] = pi_new[k1];
    }

    /*******The E-step of the EM********/

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K;  k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j + 1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    deni = pow(eps + dnorm(REAL(myY)[i] - mui, 0, sigma_new[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_new[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    phi[i][k1] = pi_new[k1] * deni;
    sumi += phi[i][k1];
    }
    for(k1=0; k1 < K;  k1++){
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    /**********End of the E-step*******/
    /**********End of the E-step*******/
    for(k1=0; k1 < K; k1++){
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    if(asInteger(myNorm) == 1){
    REAL(predict)[k1*N+i] = mui;
    REAL(residual)[k1*N+i] = REAL(myY)[i] - mui;
    }else{
    REAL(predict)[k1*N+i] = exp(mui);
    REAL(residual)[k1*N+i] = exp(REAL(myY)[i]) - exp(mui);
    }
    REAL(tau)[k1*N+i] = W[i][k1];
    }
    }

    const char *names[] = {"alpha", "beta", "sigma", "pi", "LogLik", "BIC",
       "AIC", "MaxIter", "tau", "predict", "residual",
       "selection", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, alpha);
    SET_VECTOR_ELT(res, 1, beta);
    SET_VECTOR_ELT(res, 2, sigma);
    SET_VECTOR_ELT(res, 3, pi);
    SET_VECTOR_ELT(res, 4, ScalarReal(loglike1));
    SET_VECTOR_ELT(res, 5, ScalarReal(BIC));
    SET_VECTOR_ELT(res, 6, ScalarReal(AIC));
    SET_VECTOR_ELT(res, 7, ScalarInteger(niter1));
    SET_VECTOR_ELT(res, 8, tau);
    SET_VECTOR_ELT(res, 9, predict);
    SET_VECTOR_ELT(res, 10, residual);
    SET_VECTOR_ELT(res, 11, SeL);
    UNPROTECT(9);
    return res;
}

/* ****************** Tuning Parameter Normal and Log-Normal ************** */
SEXP FMR_Norm_CTun(SEXP myY,
       SEXP myX,
       SEXP myK,
       SEXP myD,
       SEXP myN,
       SEXP mydelta,
       SEXP myAlpha,
       SEXP myBeta,
       SEXP mySigma,
       SEXP myPi,
       SEXP myacs,
       SEXP myridge,
       SEXP myEMiter,
       SEXP myeps,
       SEXP myNorm,
       SEXP myPenalty,
       SEXP myPiPor,
       SEXP myMcpTun,
       SEXP mySicaTun,
       SEXP myTol,
       SEXP myLambMin,
       SEXP myLambMax,
       SEXP mynLamb
)
{
    int i;
    int j;
    int k1;
    int l;
    int l1;
    int check1[1];

    double McpTun = asReal(myMcpTun);
    double SicaTun = asReal(mySicaTun);

    double LambMax = asReal(myLambMax);
    double LambMin = asReal(myLambMin);
    int MaxLim = asInteger(mynLamb);

    double Tol = asReal(myTol);

    int N = asInteger(myN);
    int D = asInteger(myD);
    int K = asInteger(myK);
    int jar = asInteger(myPenalty);

    int OD[K];
    for(k1=0; k1<K; k1++){
    OD[k1]=0;
    for(j=0; j<D+1; j++){
    OD[k1] += accessAcsArr(myacs, j, k1, D);       }
    }

    double ridge = asReal(myridge);
    double gam = 0.0;
    if(ridge != 0.0)
    gam = 1;

    double eps = asReal(myeps);
    double GamMP = asReal(myPiPor);

    double one_X[N][D+1][K];

    double beta_ini[D][K];
    double beta_old[D][K];
    double beta_new[D][K];

    double alpha_ini[K];
    double alpha_old[K];
    double alpha_new[K];

    double sigma_ini[K];
    double sigma_old[K];

    double pi_ini[K];

    double W[N][K];
    double phi[N][K];

    double w_s[N][K];
    double Aw[N][K];
    double V[N][K];

    double eps1[K];

    double vecder[D];
    double vecSig[D+1];
    double En[K];

    double oneXTWY[D+1];
    double oneXTWX[D+1][D+1];
    double selection[D][K];

    double oneComMat[D+2][D+2];
    double OneSolv[D+2];
    double oneComMatVec[(D+2) * (D+2)];

    double sumi;
    double mui;
    double deni;

    double n1[K];
    double BIC[MaxLim][K];
    double Max_BIC[K];
    double lambda1[MaxLim];

    int count1[K][MaxLim];
    int indx1[K];

    double hl = (LambMax - LambMin)/MaxLim;
    for(l1=0; l1 < MaxLim; l1++)
    lambda1[l1] = LambMin + l1 * (hl);

    double loglike1;

    double holdveccov[D];

    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++){
    beta_new[j][k1] = beta_old[j][k1] = beta_ini[j][k1] =
    REAL(myBeta)[D*k1+j] * accessAcsArr(myacs, j+1, k1, D);
    }
    alpha_new[k1] = alpha_old[k1] = alpha_ini[k1] =
    REAL(myAlpha)[k1] * accessAcsArr(myacs, 0, k1, D);
    /* sigma_new[k1] = */ sigma_old[k1] = sigma_ini[k1] =
    REAL(mySigma)[k1];
    /*pi_new[k1] =  pi_old[k1]  = */ pi_ini[k1] = REAL(myPi)[k1];
    }


    for(k1=0; k1 < K; k1++){
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1){
    for (i=0; i < N; i++){
    one_X[i][l][k1] = 1.0;
    }
    l++;
    }

    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j + 1, k1, D) == 1){
    for (i=0; i < N; i++){
    one_X[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    for (j=0; j < D; j++) {
    holdveccov[j] = beta_ini[j][k1];
    }
    if((jar == 1) || (jar == 2))
    eps1[k1] = eps * minimum(holdveccov, D) / (2 * N * 0.1);
    else
    eps1[k1] = eps * minimum(holdveccov, D) / (4 * N * 0.1);
    }

    for(k1=0; k1 < K; k1++)
    n1[k1] = 0.0;

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui +=  REAL(myX)[N*j+i] * beta_ini[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_ini[k1] * accessAcsArr(myacs, 0, k1, D);
    deni = pow(eps + dnorm(REAL(myY)[i] - mui, 0, sigma_ini[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_ini[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    phi[i][k1] = pi_ini[k1] * deni;
    sumi += phi[i][k1];
    }
    for(k1=0; k1 < K; k1++){
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    for(k1=0; k1 < K; k1++){
    En[k1] = N * pow(pi_ini[k1], GamMP);

    for (j=0; j < D; j++) {
    holdveccov[j] = beta_old[j][k1];
    }

    for(l1=0; l1 < MaxLim; l1++){
    if(jar == 1)
    for(j=0; j < D; j++)
    vecder[j] = lambda1[l1];
    else if (jar == 2)
    scadpen(lambda1[l1], vecder, holdveccov, D);
    else if (jar == 3)
    mcppen(lambda1[l1], vecder, holdveccov, D, McpTun);
    else if (jar == 4)
    sicapen(lambda1[l1], vecder, holdveccov, D, SicaTun);
    else if (jar == 5){
    for(j=0; j < D; j++)
    vecder[j] = lambda1[l1] / (fabs(beta_ini[j][k1]) + eps);
    }
    else
    hardpen(lambda1[l1], vecder, holdveccov, D);

    l=0;
    if(accessAcsArr(myacs, 0, k1, D)==1)
    vecSig[l++] = 0.0;
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j + 1, k1, D) == 1)
    vecSig[l++] = vecder[j] / (fabs(beta_old[j][k1]) + eps1[k1]);
    }

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < OD[k1]; j++){
    sumi = 0.0;
    for(l=0; l < N; l++)
    sumi += one_X[l][i][k1] * W[l][k1] * one_X[l][j][k1];
    if(i != j)
    oneXTWX[i][j] = sumi;
    else
    oneXTWX[i][j] = sumi + En[k1] * vecSig[j] + ridge *
    log(N);
    }
    }

    if(accessAcsArr(myacs, 0, k1, D) == 1)
    oneXTWX[0][0] = oneXTWX[0][0] - (En[k1] * vecSig[0]) - ridge * log(N);

    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    w_s[i][k1] = (REAL(myY)[i] - mui) / sigma_old[k1];
    Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 -
    pnorm(w_s[i][k1], 0, 1, 1, 0)));
    V[i][k1] = INTEGER(mydelta)[i] * REAL(myY)[i] +
    (1 - INTEGER(mydelta)[i]) * (mui + sigma_old[k1] * Aw[i][k1]);
    }

    /****** Constructing the weigthed vector oneXTWY ***/

    for(j=0; j < OD[k1]; j++){
    sumi = 0.0;
    for(i=0; i < N; i++)
    sumi += one_X[i][j][k1] * W[i][k1] * V[i][k1];
    oneXTWY[j] = sumi;
    }

    /***In a system Ax=b, adding b to A as its last column**/

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < (OD[k1] + 1); j++){
    if(j != OD[k1])
    oneComMat[i][j] = oneXTWX[i][j];
    else
    oneComMat[i][j] = oneXTWY[i];
    }
    }

    for (i=0; i < (OD[k1] + 1); i++) {
    for (j=0; j < (OD[k1] + 1); j++) {
    oneComMatVec[j+i*(OD[k1]+1)] = oneComMat[i][j];
    }
    }

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/
    /**************************************************************/
    count1[k1][l1]=0;
    sol(OD[k1], oneComMatVec, OneSolv, check1);

    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1)
    alpha_new[k1] = OneSolv[l++];
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j + 1, k1, D) == 1){
    beta_new[j][k1]  = OneSolv[l++];
    if(fabs(beta_new[j][k1]) < Tol)
    selection[j][k1]=0;
    else
    selection[j][k1] = 1;
    count1[k1][l1] += selection[j][k1];
    }
    }

    for(j=0; j < D; j++){
    beta_old[j][k1] = beta_new[j][k1];
    }

    alpha_old[k1] = alpha_new[k1];

    loglike1 = 0.0;
    n1[k1] = 0.0;

    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    selection[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D) ;
    deni = pow(eps + dnorm(REAL(myY)[i] - mui, 0, sigma_old[k1], 0),
       INTEGER(mydelta)[i]) * pow(eps + 1 - pnorm(REAL(myY)[i] - mui, 0,
       sigma_old[k1], 1, 0), 1 - INTEGER(mydelta)[i]);
    loglike1 += W[i][k1] * log(deni);
    n1[k1] += W[i][k1];
    }

    BIC[l1][k1] = loglike1 - 0.5 * (count1[k1][l1]) * log(n1[k1]) -
    gam * (count1[k1][l1]) * log(OD[k1]);

    if(l1 == 0){
    Max_BIC[k1] = BIC[l1][k1];
    indx1[k1] = l1;
    }
    else if(BIC[l1][k1] > Max_BIC[k1])
    {
    Max_BIC[k1] = BIC[l1][k1];
    indx1[k1] = l1;
    }
    }
    }

    SEXP OptLam = PROTECT(allocVector(REALSXP, K));
    for(k1=0; k1 < K; k1++){
    REAL(OptLam)[k1] = lambda1[indx1[k1]];
    }

    const char *names[] = {"OptLam", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, OptLam);
    UNPROTECT(2);
    return res;
}

/* *************************** MLE Weibull *************************** */
SEXP FMR_Weibl_MLE(SEXP myY,
       SEXP myX,
       SEXP myK,
       SEXP myD,
       SEXP myN,
       SEXP mydelta,
       SEXP myAlpha,
       SEXP myBeta,
       SEXP mySigma,
       SEXP myPi,
       SEXP myacs,
       SEXP myridge,
       SEXP myEMiter,
       SEXP myeps,
       SEXP myepsC,
       SEXP myNRiter,
       SEXP myNRPor
)
{
    int i;
    int j;
    int j1;
    int k1;
    int l;
    int niter1;
    int niter2;
    int check1[1];

    int N = asInteger(myN);
    int D = asInteger(myD);
    int K = asInteger(myK);
    int EM_Miter = asInteger(myEMiter);
    int NR_maxiter = asInteger(myNRiter);

    int OD[K];
    for(k1=0; k1<K; k1++){
    OD[k1]=0;
    for(j=0; j<D+1; j++){
    OD[k1] += accessAcsArr(myacs, j, k1, D);
    }
    }

    double eps_conv = asReal(myepsC);
    double eps = asReal(myeps);

    double ridge = asReal(myridge);
    int NRportion = asReal(myNRPor);
    int alp=0;

    double multX[N][D][K];
    double OneX[N][D+1][K];
    double sumi;
    double mui;
    double deni;

    double W[N][K];
    double phi[N][K];

    double beta_ini[D][K];
    double beta_old[D][K];
    double beta_new[D][K];

    double alpha_ini[K];
    double alpha_old[K];
    double alpha_new[K];

    double sigma_ini[K];
    double sigma_old[K];
    double sigma_new[K];

    double pi_ini[K];
    double pi_old[K];
    double pi_new[K];

    double oneXTWY[D+1];
    double oneXTWX[D+1][D+1];

    double sumwi[K];
    double sumi3[K];
    double sumi5[K];

    double SumDif;

    double OneMat[D+2][D+2];
    double OneSolv[D+2];
    double OneMatVec[(D+2) * (D+2)];
    double loglike1;
    double oldloglike1;
    double newloglike1;

    double zZ = 0.0;
    double eZ = 0.0;

    char convg1 = 'n';
    char convg2 = 'n';

    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++){
    beta_new[j][k1] = beta_old[j][k1] = beta_ini[j][k1] =
    REAL(myBeta)[D*k1+j] * accessAcsArr(myacs, j+1, k1, D);
    }
    alpha_new[k1] = alpha_old[k1] = alpha_ini[k1] =
    REAL(myAlpha)[k1] * accessAcsArr(myacs, 0, k1, D);
    sigma_new[k1] = sigma_old[k1] = sigma_ini[k1] =
    REAL(mySigma)[k1];
    pi_new[k1] = pi_old[k1] = pi_ini[k1] = REAL(myPi)[k1];
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    for (i=0; i < N; i++){
    multX[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1){
    for (i=0; i < N; i++){
    OneX[i][l][k1] = 1.0;
    }
    l++;
    }

    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    for (i=0; i < N; i++){
    OneX[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    niter1=0;
    convg1 = 'n';
    while((convg1 != 'y') && (niter1 < EM_Miter)){
    /****Beginning of each iteration****/
    /*******The E-step of the EM********/
    for (k1=0; k1 < K; k1++) {
    sumwi[k1] = 0.0;
    }

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    phi[i][k1] = pi_old[k1] * deni;
    sumi += phi[i][k1];
    }
    for(k1=0; k1 < K; k1++){
    W[i][k1] = phi[i][k1] / sumi;
    sumwi[k1] += W[i][k1];
    }
    }

    for(k1=0; k1 < K; k1++){
    pi_new[k1] = sumwi[k1] / N;
    }

    /**********End of the E-step*******/

    for (k1=0; k1 < K; k1++) {
    for (j1=0; j1 < D; j1++) {
    beta_ini[j1][k1] = beta_old[j1][k1];
    }
    alpha_ini[k1] = alpha_old[k1];
    sigma_ini[k1] = sigma_old[k1];
    pi_ini[k1] = pi_old[k1];
    }

    /*********M-step of the EM*********/

    for(k1=0; k1 < K; k1++){
    niter2=0;
    convg2 = 'n';
    while((convg2 != 'y') && (niter2 < NR_maxiter)){
    /******Compute Likelihood***/
    oldloglike1 = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) *
    exp(-eZ);
    oldloglike1 += W[i][k1] * log(deni);
    }

    /*****Constructing the weighted oneXTWY matrix***/
    j1=0;
    for(j=0; j < (D+1); j++){
    if(accessAcsArr(myacs, j, k1, D) == 1){
    sumi = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(l=0; l < D; l++)
    mui += REAL(myX)[N*l+i] * beta_old[l][k1] *
    accessAcsArr(myacs, l+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    if(j==0){
    sumi += W[i][k1] * (accessAcsArr(myacs, 0, k1, D) /
    sigma_old[k1]) * (eZ - INTEGER(mydelta)[i]);
    }
    else{
    sumi += W[i][k1] * (REAL(myX)[N*(j-1)+i] *
    accessAcsArr(myacs, j, k1, D) / sigma_old[k1]) *
    (eZ - INTEGER(mydelta)[i]);
    }
    }
    oneXTWY[j1++] =  sumi;
    }
    }

    /*****Constructing the weighted Hessian matrix***/

    for(i=0; i < OD[k1]; i++){
    for(j = i; j < OD[k1]; j++){
    sumi = 0.0;
    for(l=0; l < N; l++){
    mui = 0.0;
    for(j1=0; j1 < D; j1++)
    mui += REAL(myX)[N*j1+l] * beta_old[j1][k1] *
    accessAcsArr(myacs, j1+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[l] - mui) / sigma_old[k1]);
    sumi += - W[l][k1] * OneX[l][i][k1] * OneX[l][j][k1] /
    (sigma_old[k1] * sigma_old[k1]) * eZ;
    }
    if(i == j)
    oneXTWX[i][j] = sumi + ridge * log(N);
    else
    oneXTWX[j][i] = oneXTWX[i][j] = sumi;
    }
    }

    if(accessAcsArr(myacs, 0, k1, D) == 1)
    oneXTWX[0][0] = oneXTWX[0][0] - ridge * log(N);

    /***In a system Ax=b, adding b to A as its last column**/

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < (OD[k1] + 1); j++){
    if(j != OD[k1])
    OneMat[i][j] = - oneXTWX[i][j];
    else
    OneMat[i][j] = oneXTWY[i];
    }
    }

    for (i=0; i < (OD[k1] + 1); i++) {
    for (j=0; j < (OD[k1] + 1); j++) {
    OneMatVec[j+i*(OD[k1]+1)] = OneMat[i][j];
    }
    }

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/
    /**************************************************************/
    sol((OD[k1]), OneMatVec, OneSolv, check1);
    alp=0;
    do {
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1)
    alpha_new[k1] = pow(0.5, alp) * OneSolv[l++] +
    alpha_old[k1];
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    beta_new[j][k1] = pow(0.5, alp) *
    OneSolv[l++] + beta_old[j][k1];
    }
    }

    sumi3[k1] = 0.0;
    sumi5[k1] = 0.0;

    for(i=0; i < N; i++){
    mui = 0.0;
    for(l=0; l < D; l++)
    mui += REAL(myX)[N*l+i] * beta_new[l][k1] *
    accessAcsArr(myacs, l+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    zZ = (REAL(myY)[i] - mui) / sigma_old[k1];
    eZ = exp(zZ);
    sumi3[k1] += W[i][k1] * (- INTEGER(mydelta)[i] / sigma_old[k1] +
    (zZ / (sigma_old[k1])) * (eZ - INTEGER(mydelta)[i]));
    sumi5[k1] += W[i][k1] * (INTEGER(mydelta)[i] / (sigma_old[k1] *
    sigma_old[k1]) + (zZ / (sigma_old[k1] * sigma_old[k1])) *
    (2 * INTEGER(mydelta)[i] - (2 + zZ) * eZ));
    }
    sigma_new[k1] = sigma_old[k1] - pow(0.5, alp) * (sumi3[k1] /
    sumi5[k1]);

    newloglike1 = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_new[k1]), INTEGER(mydelta)[i]) *
    exp(-eZ);
    newloglike1 += W[i][k1] * log(deni);
    }
    alp++;

    } while ((oldloglike1 > newloglike1) & (alp < NRportion));

    SumDif = 0.0;
    niter2++;
    for(j=0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_old[j][k1], 2) *
    accessAcsArr(myacs, j+1, k1, D);
    SumDif += pow(alpha_new[k1] - alpha_old[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);
    SumDif += pow(sigma_new[k1] - sigma_old[k1], 2);
    SumDif += pow(pi_new[k1] - pi_old[k1], 2);

    convg2 = 'n';
    if(SumDif < eps_conv)
    convg2 = 'y';

    for(j=0; j < D; j++)
    beta_old[j][k1] = beta_new[j][k1];
    alpha_old[k1] = alpha_new[k1];
    sigma_old[k1] = sigma_new[k1];
    pi_old[k1] = pi_new[k1];
    }
    }

    /*****End of the M-step of the EM*****/

    SumDif = 0.0;
    niter1++;
    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_ini[j][k1], 2) *
    accessAcsArr(myacs, j+1, k1, D);
    SumDif += pow(alpha_new[k1] - alpha_ini[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);
    SumDif += pow(sigma_new[k1] - sigma_ini[k1], 2);
    SumDif += pow(pi_new[k1] - pi_ini[k1], 2);
    }

    convg1 = 'n';
    if(SumDif < eps_conv)
    convg1 = 'y';
    for(k1=0; k1 < K; k1++){
    alpha_old[k1] = alpha_new[k1];
    for(j=0; j < D; j++){
    beta_old[j][k1] = beta_new[j][k1];
    }
    sigma_old[k1] = sigma_new[k1];
    pi_old[k1] = pi_new[k1];
    }
    }

    /*******End of each iteration *******/

    loglike1 = 0.0;
    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_new[k1]);
    deni = eps + pow((eZ / sigma_new[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    phi[i][k1] = pi_new[k1] * deni;
    sumi += phi[i][k1];
    }
    loglike1 += log(sumi);
    for(k1=0; k1 < K;  k1++){
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    int SumPar = 0;
    for(k1=0; k1<K; k1++){
    SumPar += OD[k1];
    }
    SumPar += 2*K-1;

    double BIC = 0.0;
    double AIC = 0.0;
    BIC = - (2 * loglike1) + SumPar * log(N);
    AIC = - (2 * loglike1) + 2 * SumPar;

    SEXP alpha = PROTECT(allocVector(REALSXP, K));
    SEXP beta = PROTECT(allocVector(REALSXP, K*D));
    SEXP sigma = PROTECT(allocVector(REALSXP, K));
    SEXP pi = PROTECT(allocVector(REALSXP, K));
    SEXP predict = PROTECT(allocVector(REALSXP, K*N));
    SEXP residual = PROTECT(allocVector(REALSXP, K*N));
    SEXP tau = PROTECT(allocVector(REALSXP, K*N));

    for (k1=0; k1 < K; k1++){
    for (j=0;  j < D;  j++)
    REAL(beta)[k1 * D + j] = beta_new[j][k1];
    REAL(alpha)[k1] = alpha_new[k1];
    REAL(sigma)[k1] = sigma_new[k1];
    REAL(pi)[k1] = pi_new[k1];
    }

    /**********End of the E-step*******/
    for(k1=0; k1 < K; k1++){
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    REAL(predict)[k1*N+i] = exp(mui);
    REAL(residual)[k1*N+i] = exp(REAL(myY)[i]) - exp(mui);
    REAL(tau)[k1*N+i] = W[i][k1];
    }
    }
    const char *names[] = {"alpha", "beta", "sigma", "pi", "LogLik", "BIC",
       "AIC", "MaxIter", "tau", "predict", "residual", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, alpha);
    SET_VECTOR_ELT(res, 1, beta);
    SET_VECTOR_ELT(res, 2, sigma);
    SET_VECTOR_ELT(res, 3, pi);
    SET_VECTOR_ELT(res, 4, ScalarReal(loglike1));
    SET_VECTOR_ELT(res, 5, ScalarReal(BIC));
    SET_VECTOR_ELT(res, 6, ScalarReal(AIC));
    SET_VECTOR_ELT(res, 7, ScalarInteger(niter1));
    SET_VECTOR_ELT(res, 8, tau);
    SET_VECTOR_ELT(res, 9, predict);
    SET_VECTOR_ELT(res, 10, residual);
    UNPROTECT(8);
    return res;
}

/* ************************* Variable Selection Weibull***************** */
SEXP FMR_Weibl_MPLE(SEXP myY,
    SEXP myX,
    SEXP myK,
    SEXP myD,
    SEXP myN,
    SEXP mydelta,
    SEXP myAlpha,
    SEXP myBeta,
    SEXP mySigma,
    SEXP myPi,
    SEXP myacs,
    SEXP myridge,
    SEXP myEMiter,
    SEXP myeps,
    SEXP myepsC,
    SEXP myNRiter,
    SEXP myPenalty,
    SEXP mylambda,
    SEXP myPiPor,
    SEXP myMcpTun,
    SEXP mySicaTun,
    SEXP myTol,
    SEXP myNRPor
)
{
    int i;
    int j;
    int j1;
    int k1;
    int l;
    int niter1;
    int niter2;
    int check1[1];

    int N = asInteger(myN);
    int D = asInteger(myD);
    int K = asInteger(myK);
    int jar = asInteger(myPenalty);
    double McpTun = asReal(myMcpTun);
    double SicaTun = asReal(mySicaTun);

    double Tol = asReal(myTol);

    int EM_Miter = asInteger(myEMiter);
    int NR_maxiter = asInteger(myNRiter);

    int OD[K];
    for(k1=0; k1<K; k1++){
    OD[k1]=0;
    for(i=0; i<D+1; i++){
    OD[k1] += accessAcsArr(myacs, i, k1, D);
    }
    }

    double eps_conv = asReal(myepsC);
    double eps = asReal(myeps);
    double GamMP = asReal(myPiPor);

    double ridge = asReal(myridge);
    int NRportion = asReal(myNRPor);
    int alp=0;

    double multX[N][D][K];
    double OneX[N][D+1][K];
    double sumi;
    double mui;
    double deni;

    double W[N][K];
    double phi[N][K];

    double optlam[K];

    double beta_ini[D][K];
    double beta_old[D][K];
    double beta_new[D][K];

    double alpha_ini[K];
    double alpha_old[K];
    double alpha_new[K];

    double sigma_ini[K];
    double sigma_old[K];
    double sigma_new[K];

    double pi_ini[K];
    double pi_old[K];
    double pi_new[K];

    double oneXTWY[D+1];
    double oneXTWX[D+1][D+1];
    double select[D][K];

    double eps1[K];
    double vecder[D];
    double vecF[D+1];
    double En[K];

    double sumwi[K];
    double sumi3[K];
    double sumi5[K];

    double SumDif;

    double OneMat[D+2][D+2];
    double OneSolv[D+2];
    double OneMatVec[(D+2) * (D+2)];
    double loglike1;
    double oldloglike1;
    double newloglike1;

    double Hcov[D];

    for (k1=0; k1 < K; k1++) {
    optlam[k1]= REAL(mylambda)[k1];
    }

    double zZ = 0.0;
    double eZ = 0.0;

    char convg1 = 'n';
    char convg2 = 'n';

    for(k1=0; k1 < K; k1++){
    for(i=0; i < D; i++){
    beta_new[i][k1] = beta_old[i][k1] =
    beta_ini[i][k1] = REAL(myBeta)[D*k1+i] * accessAcsArr(myacs, i+1, k1, D);
    }
    alpha_new[k1] = alpha_old[k1] = alpha_ini[k1] =
    REAL(myAlpha)[k1] * accessAcsArr(myacs, 0, k1, D);
    sigma_new[k1] = sigma_old[k1] = sigma_ini[k1] =
    REAL(mySigma)[k1];
    pi_new[k1] = pi_old[k1] = pi_ini[k1] = REAL(myPi)[k1];
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    for (i=0; i < N; i++){
    multX[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1){
    for (i=0; i < N; i++){
    OneX[i][l][k1] = 1.0;
    }
    l++;
    }

    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    for (i=0; i < N; i++){
    OneX[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    for (i=0; i < D; i++) {
    Hcov[i] = beta_ini[i][k1];
    }
    if((jar == 1) || (jar == 2))
    eps1[k1] = eps * minimum(Hcov, D) / (2 * N * optlam[k1]);
    else
    eps1[k1] = eps * minimum(Hcov, D) / (4 * N * optlam[k1]);
    }

    niter1=0;
    convg1 = 'n';
    while((convg1 != 'y') && (niter1 < EM_Miter)){
    /****Beginning of each iteration****/
    /*******The E-step of the EM********/
    for (k1=0; k1 < K; k1++) {
    sumwi[k1] = 0.0;
    }

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    phi[i][k1] = pi_old[k1] * deni;
    sumi += phi[i][k1];
    }
    for(k1=0; k1 < K; k1++){
    W[i][k1] = phi[i][k1] / sumi;
    sumwi[k1] += W[i][k1];
    }
    }

    for(k1=0; k1 < K; k1++){
    pi_new[k1] = sumwi[k1] / N;
    }

    /**********End of the E-step*******/

    for (k1=0; k1 < K; k1++) {
    for (j1=0; j1 < D; j1++) {
    beta_ini[j1][k1] = beta_old[j1][k1];
    }
    alpha_ini[k1] = alpha_old[k1];
    sigma_ini[k1] = sigma_old[k1];
    pi_ini[k1] = pi_old[k1];
    }

    /*********M-step of the EM*********/

    for(k1=0; k1 < K; k1++){

    En[k1] = N * pow(pi_old[k1], GamMP);

    for (j=0; j < D; j++) {
    Hcov[j] = beta_old[j][k1];
    }

    if(jar == 1)
    for(j=0; j < D; j++)
    vecder[j] = optlam[k1];
    else if (jar == 2)
    scadpen(optlam[k1], vecder, Hcov, D);
    else if (jar == 3)
    mcppen(optlam[k1], vecder, Hcov, D, McpTun);
    else if (jar == 4)
    sicapen(optlam[k1], vecder, Hcov, D, SicaTun);
    else if (jar == 5){
    for(j=0; j < D; j++)
    vecder[j] = optlam[k1] / (fabs(beta_ini[j][k1]) + eps);
    }
    else
    hardpen(optlam[k1], vecder, Hcov, D);


    l=0;
    if(accessAcsArr(myacs, l, k1, D) == 1)
    vecF[l++] = 0.0;
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1)
    vecF[l++] = vecder[j] / (fabs(beta_old[j][k1]) + eps1[k1]);
    }

    niter2=0;
    convg2 = 'n';
    while((convg2 != 'y') && (niter2 < NR_maxiter)){
    /******Compute Likelihood***/
    oldloglike1 = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    oldloglike1 += W[i][k1] * log(deni);
    }

    /*****Constructing the weighted oneXTWY matrix***/

    j1=0;
    for(j=0; j < (D+1); j++){
    if(accessAcsArr(myacs, j, k1, D) == 1){
    sumi = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(l=0; l < D; l++)
    mui += REAL(myX)[N*l+i] * beta_old[l][k1] * accessAcsArr(myacs, l+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    if(j==0){
    sumi += W[i][k1] * (accessAcsArr(myacs, 0, k1, D) / sigma_old[k1]) *
    (eZ - INTEGER(mydelta)[i]);
    }
    else{
    sumi += W[i][k1] * (REAL(myX)[N*(j-1)+i] *
    accessAcsArr(myacs, j, k1, D) /
    sigma_old[k1]) * (eZ - INTEGER(mydelta)[i]);
    }
    }
    oneXTWY[j1++] =  sumi;
    }
    }

    /*****Constructing the weighted Hessian matrix***/

    for(i=0; i < OD[k1]; i++){
    for(j = i; j < OD[k1]; j++){
    sumi = 0.0;
    for(l=0; l < N; l++){
    mui = 0.0;
    for(j1=0; j1 < D; j1++)
    mui += REAL(myX)[N*j1+l] * beta_old[j1][k1] *
    accessAcsArr(myacs, j1+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[l] - mui) / sigma_old[k1]);
    sumi += - W[l][k1] * OneX[l][i][k1] * OneX[l][j][k1] /
    (sigma_old[k1] * sigma_old[k1]) * eZ;
    }
    if(i == j)
    oneXTWX[i][j] = sumi - En[k1] * vecF[j] + ridge * log(N);
    else
    oneXTWX[j][i] = oneXTWX[i][j] = sumi;
    }
    }

    if(accessAcsArr(myacs, 0, k1, D) == 1)
    oneXTWX[0][0] = oneXTWX[0][0] + En[k1]* vecF[0] - ridge * log(N);

    /***In a system Ax=b, adding b to A as its last column**/

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < (OD[k1] + 1); j++){
    if(j != OD[k1])
    OneMat[i][j] = - oneXTWX[i][j];
    else
    OneMat[i][j] = oneXTWY[i];
    }
    }

    for (i=0; i < (OD[k1] + 1); i++) {
    for (j=0; j < (OD[k1] + 1); j++) {
    OneMatVec[j+i*(OD[k1]+1)] = OneMat[i][j];
    }
    }

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/
    /**************************************************************/
    sol((OD[k1]), OneMatVec, OneSolv, check1);
    alp=0;
    do {
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1)
    alpha_new[k1] = pow(0.5, alp) * OneSolv[l++] +
    alpha_old[k1];
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    beta_new[j][k1] = pow(0.5, alp) *
    OneSolv[l++] + beta_old[j][k1];
    if(fabs(beta_new[j][k1]) <= Tol){
    beta_new[j][k1] = 0.0;
    }

    }
    }

    sumi3[k1] = 0.0;
    sumi5[k1] = 0.0;

    for(i=0; i < N; i++){
    mui = 0.0;
    for(l=0; l < D; l++)
    mui += REAL(myX)[N*l+i] * beta_new[l][k1] *
    accessAcsArr(myacs, l+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    zZ = (REAL(myY)[i] - mui) / sigma_old[k1];
    eZ = exp(zZ);
    sumi3[k1] += W[i][k1] * (- INTEGER(mydelta)[i] / sigma_old[k1] +
    (zZ / (sigma_old[k1])) * (eZ - INTEGER(mydelta)[i]));
    sumi5[k1] += W[i][k1] * (INTEGER(mydelta)[i] / (sigma_old[k1] *
    sigma_old[k1]) + (zZ / (sigma_old[k1] * sigma_old[k1])) *
    (2 * INTEGER(mydelta)[i] - (2 + zZ) * eZ));
    }
    sigma_new[k1] = sigma_old[k1] - pow(0.5, alp) * (sumi3[k1] /
    sumi5[k1]);

    newloglike1 = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_new[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    newloglike1 += W[i][k1] * log(deni);
    }
    alp++;

    } while ((oldloglike1 > newloglike1) & (alp < NRportion));

    SumDif = 0.0;
    niter2++;
    for(j=0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_old[j][k1], 2) *
    accessAcsArr(myacs, j+1, k1, D);
    SumDif += pow(alpha_new[k1] - alpha_old[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);
    SumDif += pow(sigma_new[k1] - sigma_old[k1], 2);
    SumDif += pow(pi_new[k1] - pi_old[k1], 2);

    convg2 = 'n';
    if(SumDif < eps_conv)
    convg2 = 'y';

    for(j=0; j < D; j++)
    beta_old[j][k1] = beta_new[j][k1];
    alpha_old[k1] = alpha_new[k1];
    sigma_old[k1] = sigma_new[k1];
    pi_old[k1] = pi_new[k1];
    }
    }

    /*****End of the M-step of the EM*****/

    SumDif = 0.0;
    niter1++;
    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_ini[j][k1], 2) *
    accessAcsArr(myacs, j+1, k1, D);
    SumDif += pow(alpha_new[k1] - alpha_ini[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);
    SumDif += pow(sigma_new[k1] - sigma_ini[k1], 2);
    SumDif += pow(pi_new[k1] - pi_ini[k1], 2);
    }

    convg1 = 'n';
    if(SumDif < eps_conv)
    convg1 = 'y';
    for(k1=0; k1 < K; k1++){
    alpha_old[k1] = alpha_new[k1];
    for(j=0; j < D; j++){
    beta_old[j][k1] = beta_new[j][k1];
    }
    sigma_old[k1] = sigma_new[k1];
    pi_old[k1] = pi_new[k1];
    }
    }

    /*******End of each iteration *******/

    loglike1 = 0.0;
    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_new[k1]);
    deni = eps + pow((eZ / sigma_new[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    phi[i][k1] = pi_new[k1] * deni;
    sumi += phi[i][k1];
    }
    loglike1 += log(sumi);
    for(k1=0; k1 < K;  k1++){
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    //*********************************************************************
    //*Storing the estimates of the regression coefficients
    //*in a global variable called "Betahat".
    //*Selecting the finial model and storing
    //*********************************************************************

    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++)
    select[j][k1] = 0;
    }

    for(k1=0; k1 < K; k1++){
    for(j=0; j < D; j++){
    if(fabs(beta_new[j][k1]) <= Tol){
    beta_new[j][k1] = 0.0;
    select[j][k1]= 0;
    }
    else
    select[j][k1] = 1;
    }
    }

    loglike1 = 0.0;
    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_new[k1]);
    deni = eps + pow((eZ / sigma_new[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    phi[i][k1] = pi_new[k1] * deni;
    sumi += phi[i][k1];
    }
    loglike1 += log(sumi);
    for(k1=0; k1 < K;  k1++){
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    int SUM1 = 0;
    for(k1=0; k1 < K;  k1++){
    SUM1 += accessAcsArr(myacs, 0, k1, D);
    for(j=0; j < D; j++){
    SUM1 += select[j][k1];
    }
    }

    SUM1 += 2*K-1;
    double BIC = 0.0;
    double AIC = 0.0;
    BIC = -(2 * loglike1) + SUM1 * log(N);
    AIC = -(2 * loglike1) + 2* SUM1;

    SEXP alpha = PROTECT(allocVector(REALSXP, K));
    SEXP beta = PROTECT(allocVector(REALSXP, K*D));
    SEXP sigma = PROTECT(allocVector(REALSXP, K));
    SEXP pi = PROTECT(allocVector(REALSXP, K));
    SEXP predict = PROTECT(allocVector(REALSXP, K*N));
    SEXP residual = PROTECT(allocVector(REALSXP, K*N));
    SEXP tau = PROTECT(allocVector(REALSXP, K*N));
    SEXP SeL = PROTECT(allocVector(INTSXP, K*D));

    for (k1=0; k1 < K; k1++){
    for (j=0;  j < D;  j++)
    INTEGER(SeL)[k1 * D + j] = select[j][k1];
    }

    for (k1=0; k1 < K; k1++){
    for (j=0;  j < D;  j++)
    {
    REAL(beta)[k1 * D + j] = beta_new[j][k1];
    }
    REAL(alpha)[k1] = alpha_new[k1];
    REAL(sigma)[k1] = sigma_new[k1];
    REAL(pi)[k1] = pi_new[k1];
    }

    /**********End of the E-step*******/
    for(k1=0; k1 < K; k1++){
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += multX[i][j][k1] * beta_new[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    REAL(predict)[k1*N+i] = exp(mui);
    REAL(residual)[k1*N+i] = exp(REAL(myY)[i]) - exp(mui);
    REAL(tau)[k1*N+i] = W[i][k1];
    }
    }

    const char *names[] = {"alpha", "beta", "sigma", "pi", "LogLik", "BIC",
       "AIC", "MaxIter", "tau", "predict", "residual",
       "selection", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, alpha);
    SET_VECTOR_ELT(res, 1, beta);
    SET_VECTOR_ELT(res, 2, sigma);
    SET_VECTOR_ELT(res, 3, pi);
    SET_VECTOR_ELT(res, 4, ScalarReal(loglike1));
    SET_VECTOR_ELT(res, 5, ScalarReal(BIC));
    SET_VECTOR_ELT(res, 6, ScalarReal(AIC));
    SET_VECTOR_ELT(res, 7, ScalarInteger(niter1));
    SET_VECTOR_ELT(res, 8, tau);
    SET_VECTOR_ELT(res, 9, predict);
    SET_VECTOR_ELT(res, 10, residual);
    SET_VECTOR_ELT(res, 11, SeL);
    UNPROTECT(9);
    return res;
}

/* *************************** Tuning Parameter Weibull ******************* */
SEXP FMR_Weibl_CTun(SEXP myY,
    SEXP myX,
    SEXP myK,
    SEXP myD,
    SEXP myN,
    SEXP mydelta,
    SEXP myAlpha,
    SEXP myBeta,
    SEXP mySigma,
    SEXP myPi,
    SEXP myacs,
    SEXP myridge,
    SEXP myEMiter,
    SEXP myeps,
    SEXP myNRiter,
    SEXP myPenalty,
    SEXP myPiPor,
    SEXP myMcpTun,
    SEXP mySicaTun,
    SEXP myTol,
    SEXP myLambMin,
    SEXP myLambMax,
    SEXP mynLamb,
    SEXP myNRPor
)
{
    int i;
    int j;
    int j1;
    int k1;
    int l;
    int check1[1];
    int l1;
    int MaxLim = asInteger(mynLamb);

    double eZ = 0.0;

    int N = asInteger(myN);
    int D = asInteger(myD);
    int K = asInteger(myK);
    int jar = asInteger(myPenalty);

    double McpTun = asReal(myMcpTun);
    double SicaTun = asReal(mySicaTun);

    double Tol = asReal(myTol);

    double LambMax = asReal(myLambMax);
    double LambMin = asReal(myLambMin);

    int OD[K];
    for(k1=0; k1<K; k1++){
    OD[k1]=0;
    for(i=0; i<D+1; i++){
    OD[k1] += accessAcsArr(myacs, i, k1, D);
    }
    }

    double ridge = asReal(myridge);
    double gam = 0.0;
    if(ridge != 0.0)
    gam = 1;

    double eps = asReal(myeps);
    double GamMP = asReal(myPiPor);

    double NR_maxiter = asReal(myNRiter);
    int NRportion = asReal(myNRPor);
    int alp;

    double multX[N][D][k1];
    double OneX[N][D+1][k1];

    double optlam[K];

    double beta_ini[D][K];
    double beta_old[D][K];
    double beta_new[D][K];

    double alpha_ini[K];
    double alpha_old[K];
    double alpha_new[K];

    double sigma_ini[K];
    double sigma_old[K];
    double sigma_new[K];

    double pi_ini[K];
    double pi_old[K];
    double pi_new[K];

    double W[N][K];
    double phi[N][K];

    double eps1[K];

    double vecder[D];
    double vecF[D+1];
    double En[K];

    double oneXTWY[D+1];
    double oneXTWX[D+1][D+1];
    double select[D][K];

    double OneMat[D+2][D+2];
    double OneSolv[D+2];
    double OneMatVec[(D+2) * (D+2)];

    double sumi;
    double mui;
    double deni;

    double n1[K];
    double BIC[MaxLim][K];
    double Max_BIC[K];
    double lambda1[MaxLim];
    double newloglike1;
    double oldloglike1;
    double SumDif;

    int count1[K][MaxLim];
    int indx1[K];

    int niter2;
    char convg2;

    double hl = (LambMax - LambMin) / MaxLim;
    for(l1=0; l1 < MaxLim; l1++)
    lambda1[l1] = LambMin + l1 * (hl);

    double loglike1;

    double Hcov[D];

    for(k1=0; k1 < K; k1++){
    for(i=0; i < D; i++){
    beta_new[i][k1] = beta_old[i][k1] = beta_ini[i][k1] =
    REAL(myBeta)[D*k1+i] * accessAcsArr(myacs, i+1, k1, D);
    }
    alpha_new[k1] = alpha_old[k1] = alpha_ini[k1] =
    REAL(myAlpha)[k1] * accessAcsArr(myacs, 0, k1, D);
    sigma_new[k1] = sigma_old[k1] = sigma_ini[k1] =
    REAL(mySigma)[k1];
    pi_new[k1] = pi_old[k1] = pi_ini[k1] = REAL(myPi)[k1];
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    for (i=0; i < N; i++){
    multX[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1){
    for (i=0; i < N; i++){
    OneX[i][l][k1] = 1.0;
    }
    l++;
    }

    for (j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    for (i=0; i < N; i++){
    OneX[i][l][k1] = REAL(myX)[N*j+i];
    }
    l++;
    }
    }
    }

    for(k1=0; k1 < K; k1++){
    for (i=0; i < D; i++) {
    Hcov[i] = beta_ini[i][k1];
    }
    if((jar == 1) || (jar == 2))
    eps1[k1] = eps * minimum(Hcov, D) / (2 * N * 0.1);
    else
    eps1[k1] = eps * minimum(Hcov, D) / (4 * N * 0.1);
    }


    for(k1=0; k1 < K; k1++)
    n1[k1] = 0.0;

    for(i=0; i < N; i++){
    sumi = 0.0;
    for(k1=0; k1 < K; k1++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    phi[i][k1] = pi_old[k1] * deni;
    sumi += phi[i][k1];
    }
    for(k1=0; k1 < K; k1++){
    W[i][k1] = phi[i][k1] / sumi;
    }
    }

    for(k1=0; k1 < K; k1++){

    En[k1] = N * pow(pi_old[k1], GamMP);

    for(l1=0; l1 < MaxLim; l1++){

    niter2=0;
    convg2 = 'n';
    while((convg2 != 'y') && (niter2 < NR_maxiter)){

    oldloglike1 = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_old[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    oldloglike1 += W[i][k1] * log(deni);
    }

    for (j=0; j < D; j++) {
    Hcov[j] = beta_old[j][k1];
    }

    if(jar == 1)
    for(j=0; j < D; j++)
    vecder[j] = lambda1[l1];
    else if (jar == 2)
    scadpen(lambda1[l1], vecder, Hcov, D);
    else if (jar == 3)
    mcppen(lambda1[l1], vecder, Hcov, D, McpTun);
    else if (jar == 4)
    sicapen(lambda1[l1], vecder, Hcov, D, SicaTun);
    else if (jar == 5){
    for(j=0; j < D; j++)
    vecder[j] = optlam[k1] / (fabs(beta_ini[j][k1]) + eps);
    }
    else
    hardpen(optlam[k1], vecder, Hcov, D);

    l=0;
    if(accessAcsArr(myacs, l, k1, D)==1)
    vecF[l++] = 0.0;
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1)
    vecF[l++] = vecder[j] / (fabs(beta_old[j][k1]) +
    eps1[k1]);
    }

    /******Constructing the Hessian matrix H ***/
    for(i=0; i < OD[k1]; i++){
    for(j = i; j < OD[k1]; j++){
    sumi = 0.0;
    for(l=0; l < N; l++){
    mui = 0.0;
    for(j1=0; j1 < D; j1++)
    mui += REAL(myX)[N*j1+l] * beta_old[j1][k1] *
    accessAcsArr(myacs, j1+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[l] - mui) / sigma_old[k1]);
    sumi += - W[l][k1] * OneX[l][i][k1] * OneX[l][j][k1] /
    (sigma_old[k1] * sigma_old[k1]) * eZ;
    }
    if(i == j)
    oneXTWX[i][j] = sumi - En[k1] * vecF[j] + ridge * log(N);
    else
    oneXTWX[j][i] = oneXTWX[i][j] = sumi;
    }
    }

    if(accessAcsArr(myacs, 0, k1, D) == 1)
    oneXTWX[0][0] = oneXTWX[0][0] + En[k1] * vecF[0] - ridge * log(N);

    /******Constructing the weighted vector XTWY***/

    j1=0;
    for(j=0; j < (D+1); j++){
    if(accessAcsArr(myacs, j, k1, D) == 1){
    sumi = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(l=0; l < D; l++)
    mui += REAL(myX)[N*l+i] * beta_old[l][k1] * accessAcsArr(myacs, l+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    if(j==0){
    sumi += W[i][k1] * (accessAcsArr(myacs, 0, k1, D) / sigma_old[k1]) *
    (eZ - INTEGER(mydelta)[i]);
    }
    else{
    sumi += W[i][k1] * (REAL(myX)[N*(j-1)+i] *
    accessAcsArr(myacs, j, k1, D) /
    sigma_old[k1]) * (eZ - INTEGER(mydelta)[i]) -
    En[k1] * vecF[j1] * beta_old[j-1][k1];
    }
    }
    oneXTWY[j1++] =  sumi;
    }
    }

    /***In a system Ax=b, adding b to A as its last column**/

    for(i=0; i < OD[k1]; i++){
    for(j=0; j < (OD[k1] + 1); j++){
    if(j != OD[k1])
    OneMat[i][j] = - oneXTWX[i][j];
    else
    OneMat[i][j] = oneXTWY[i];
    }
    }

    for (i=0; i < (OD[k1] + 1); i++) {
    for (j=0; j < (OD[k1] + 1); j++) {
    OneMatVec[j+i*(OD[k1]+1)] = OneMat[i][j];
    }
    }

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/
    /**************************************************************/
    count1[k1][l1]=0;
    sol(OD[k1] , OneMatVec, OneSolv, check1);

    alp=0;
    do {
    l=0;
    if(accessAcsArr(myacs, 0, k1, D) == 1)
    alpha_new[k1] = pow(0.5, alp) * OneSolv[l++] +
    alpha_old[k1];
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    beta_new[j][k1] = pow(0.5, alp) *
    OneSolv[l++] + beta_old[j][k1];
    }
    }

    newloglike1 = 0.0;
    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_new[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_new[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    newloglike1 += W[i][k1] * log(deni);
    }
    alp++;

    } while ((oldloglike1 > newloglike1) & (alp < NRportion));

    SumDif = 0.0;
    niter2++;

    for(j=0; j < D; j++)
    SumDif += pow(beta_new[j][k1] - beta_old[j][k1], 2) *
    accessAcsArr(myacs, j+1, k1, D);
    SumDif +=pow(alpha_new[k1] - alpha_old[k1], 2) *
    accessAcsArr(myacs, 0, k1, D);

    convg2 = 'n';
    if (SumDif < eps)
    convg2 = 'y';

    count1[k1][l1]=0;
    for(j=0; j < D; j++){
    if(accessAcsArr(myacs, j+1, k1, D) == 1){
    if(fabs(beta_new[j][k1]) < Tol)
    select[j][k1]=0;
    else{
    select[j][k1] = 1;
    count1[k1][l1] += 1;
    }
    }else{
    select[j][k1]=0;
    }
    }

    for(j=0; j < D; j++)
    beta_old[j][k1] = beta_new[j][k1];
    alpha_old[k1] = alpha_new[k1];
    sigma_old[k1] = sigma_new[k1];
    }

    loglike1 = 0.0;
    n1[k1] = 0.0;

    for(i=0; i < N; i++){
    mui = 0.0;
    for(j=0; j < D; j++)
    mui += REAL(myX)[N*j+i] * beta_new[j][k1] *
    select[j][k1] * accessAcsArr(myacs, j+1, k1, D);
    mui += alpha_old[k1] * accessAcsArr(myacs, 0, k1, D);
    eZ = exp((REAL(myY)[i] - mui) / sigma_old[k1]);
    deni = eps + pow((eZ / sigma_old[k1]), INTEGER(mydelta)[i]) * exp(-eZ);
    loglike1 += W[i][k1] * log(deni);
    n1[k1] += W[i][k1];
    }

    BIC[l1][k1] = loglike1 - 0.5 * (count1[k1][l1]) * log(n1[k1]) -
    gam * (count1[k1][l1]) * log(OD[k1]);

    if(l1 == 0){
    Max_BIC[k1] = BIC[l1][k1];
    indx1[k1] = l1;
    }
    else if(BIC[l1][k1] > Max_BIC[k1])
    {
    Max_BIC[k1] = BIC[l1][k1];
    indx1[k1] = l1;
    }
    }
    }

    SEXP OptLam = PROTECT(allocVector(REALSXP, K));
    for(k1=0; k1 < K; k1++){
    REAL(OptLam)[k1] = lambda1[indx1[k1]];
    }

    const char *names[] = {"OptLam", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, OptLam);
    UNPROTECT(2);
    return res;
}
