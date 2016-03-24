/* FMRs Models; MLE, Tuning Parameter and Variable Selection */

/* ******************** Loading Liberaries ************************* */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* ******************** List of Global Variables ******************* */

#define gamma1 4.5
#define gamma2 5

/* ******************** Defining Necessary Functions *************** */

void backsub(int n, double *a, double *y)
{
  int i, k, j;
  double atemp[n + 1][n + 1];

  for (i = 0; i < (n + 1); i++) {
    for (j = 0; j < (n + 1); j++) {
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
  double max = 0;
  int loc_maxrow, loc_checksin;
  int i, j;
  loc_maxrow = row;
  loc_checksin = 1;
  double atemp[n + 1][n + 1];

  for (i = 0; i < (n + 1); i++) {
    for (j = 0; j < (n + 1); j++) {
      atemp[i][j] = a[j + i * (n + 1)];
    }
  }

  for (i = row; i < n; i++) {
    if ( fabs(atemp[i][row]) > max) {
      max = fabs(atemp[i][row]);
      loc_maxrow = i;
    }

    if (max == 0) {
      loc_checksin = 0;
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
  k = 0;

  for (i = 0; i < (n + 1); i++) {
    for (j = 0; j < (n + 1); j++) {
      atemp[i][j] = a[j + i * (n + 1)];
    }
  }

  for (i = 0; i < n; i++) {
    maxabs(i, n, a, checksin, k_vec);
    k = k_vec[0];

    if (checksin[0] == 0) {
      break;
    }

    if ( k != i)
      for (j = i; j < n + 1; j++) {
        temp = atemp[i][j];
        atemp[i][j] = atemp[k][j];
        atemp[k][j] = temp;
      } /*for j */

for (j = i + 1; j < n; j++) {
  temp = atemp[j][i] / atemp[i][i];
  for (k = i; k < n + 1; k++)
    atemp[j][k] = atemp[j][k] - 1.0 * temp * atemp[i][k];
} /* for j */
  } /* for i */

for (i = 0; i < (n + 1); i++) {
  for (j = 0; j < (n + 1); j++) {
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

double minimum(double *vector1, int ncov1)
{
  int i;
  double min1;
  min1 = fabs(vector1[0]);

  for(i = 0; i < ncov1; i++)
    if(fabs(vector1[i]) <= min1)
      min1 = fabs(vector1[i]);

    return(min1);
}
/* ************************************************************************ */

double maximum(double *vector1, int ncov1)
{
  int i;
  double max1;
  max1 = fabs(vector1[0]);

  for(i = 0; i < ncov1; i++)
    if(fabs(vector1[i]) >= max1)
      max1 = fabs(vector1[i]);

    return(max1);
}

/* ************************************************************************ */

void sicapen(double lam1, double *moshtagh, double *regcoef, int ncov1)
{
  int j;

  for(j = 0; j < ncov1; j++)
    moshtagh[j] = lam1 * gamma2 * (gamma2 + 1) / pow(gamma2 + fabs(regcoef[j]), 2);

}

/* ************************************************************************ */

void mcppen(double lam1, double *moshtagh, double *regcoef, int ncov1)
{
  int j;

  for(j = 0; j < ncov1; j++){
    if(fabs(regcoef[j]) <= (lam1 * gamma1))
      moshtagh[j] = lam1 * (1 - fabs(regcoef[j]) / (lam1 * gamma1));
    else
      moshtagh[j]= 0.0;
  }

}

/* ************************************************************************ */

void scadpen(double lam1, double *moshtagh, double *regcoef, int ncov1)
{
  int j;
  double a1;
  a1=3.7;

  for(j = 0; j < ncov1; j++){
    if(fabs(regcoef[j]) <= lam1)
      moshtagh[j] = lam1;
    else if((a1 * lam1 - fabs(regcoef[j])) <= 0)
      moshtagh[j] = 0;
    else
      moshtagh[j]= (a1 * lam1 - fabs(regcoef[j])) / (a1 - 1);
  }
}

/* ************************************************************************ */

void hardpen(double lam1, double *moshtagh, double *regcoef, int ncov1)
{
  int j;

  for(j = 0; j < ncov1; j++){
    if(fabs(regcoef[j]) <= lam1)
      moshtagh[j]= (-2) * (fabs(regcoef[j]) - lam1);
    else
      moshtagh[j] = 0;
  }
}

/* *************************** MLE Normal and Log-Normal ******************* */

#ifdef __cplusplus
extern "C" {
#endif
  void FMR_Norm_Surv_EM_MLE(double *resp,
                            double *myX,
                            double *delta,
                            double *myridgepen_prop,
                            int *myNCOMP,
                            int *myNCOV,
                            int *mynsize,
                            int *myEMmaxiter,
                            int *MaxEMiter,
                            double *myinitial_alpha,
                            double *myinitial_beta,
                            double *myinitial_sigma,
                            double *myinitial_pi,
                            double *myeps,
                            double *myepsconv,
                            double *myalpha,
                            double *mybeta,
                            double *mysigma,
                            double *mypi,
                            double *loglike,
                            double *BIC,
                            double *AIC,
                            double *GCV,
                            double *EBIC1,
                            double *EBIC5,
                            double *GIC,
                            double *predict,
                            double *residual,
                            double *tau
  )
  {
    int nsize = *mynsize;
    int NCOV = *myNCOV;
    int EM_maxiter = *myEMmaxiter;
    int NCOMP = *myNCOMP;

    double eps_conv = *myepsconv;
    double eps = *myeps;

    double ridge1 = *myridgepen_prop;

    int i;
    int j;
    int k1;
    int l;
    int niter1;
    int check1[1];

    double multX[nsize][NCOV];
    double one_X[nsize][NCOV + 1];
    double sumi;
    double mui;
    double deni;

    double W[nsize][NCOMP];
    double phi[nsize][NCOMP];

    double initbeta[NCOV][NCOMP];
    double beta0[NCOV][NCOMP];
    double new_beta0[NCOV][NCOMP];
    double beta0hat[NCOV][NCOMP];

    double initalpha[NCOMP];
    double alpha0[NCOMP];
    double new_alpha0[NCOMP];
    double alpha0hat[NCOMP];

    double initsigma[NCOMP];
    double sigma0[NCOMP];
    double new_sigma0[NCOMP];
    double sigma0hat[NCOMP];

    double sigpennom = 0.0;
    double sigpendenom = 0.0;

    double initpi[NCOMP];
    double pi0[NCOMP];
    double new_pi0[NCOMP];
    double pi0hat[NCOMP];

    double w_s[nsize][NCOMP];
    double Aw[nsize][NCOMP];
    double V[nsize][NCOMP];

    double oneXTWY[NCOV + 1];
    double oneXTWX[NCOV + 1][NCOV + 1];

    double sumwi[NCOMP];
    double sumi5[NCOMP];

    double jamconvg1;

    double oneComMat[NCOV + 2][NCOV + 2];
    double onesolution1[NCOV + 2];
    double oneComMatVec[(NCOV + 2) * (NCOV + 2)];
    double loglike1;

    char convg1 = 'n';

    double sumi1, sumi2;
    double sumi3[NCOMP];

    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < NCOV; i++){
        new_beta0[i][k1] = beta0hat[i][k1] = beta0[i][k1] = initbeta[i][k1] = myinitial_beta[NCOV * k1 + i];
      }
      new_alpha0[k1] = alpha0hat[k1] = alpha0[k1] = initalpha[k1] = myinitial_alpha[k1];
      new_sigma0[k1] = sigma0hat[k1] = sigma0[k1] = initsigma[k1] = myinitial_sigma[k1];
      new_pi0[k1] = pi0hat[k1] = pi0[k1] = initpi[k1] = myinitial_pi[k1];
    }

    for (j = 0; j < NCOV; j++){
      for (i = 0; i < nsize; i++){
        multX[i][j] = myX[nsize * j + i];
      }
    }

    for(i = 0; i < nsize; i++){
      one_X[i][0] = 1.0;
      for(j = 1; j < (NCOV + 1); j++){
        one_X[i][j] = multX[i][j - 1];
      }
    }

    niter1 = 0;
    while((convg1 != 'y') && (niter1 < EM_maxiter)){

      /****Beginning of each iteration****/

      /*******The E-step of the EM********/

      for(k1 = 0; k1 < NCOMP; k1++){
        sumwi[k1] = 0.0;
        sumi5[k1] = 0.0;
        sumi3[k1] = 0.0;
      }

      for(i = 0; i < nsize; i++){
        sumi = 0.0;
        for(k1 = 0; k1 < NCOMP;  k1++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * beta0[j][k1];
          mui += alpha0[k1];
          deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0[k1], 1, 0), 1 - delta[i]) ;
          phi[i][k1] = pi0[k1] * deni;
          sumi += phi[i][k1];
        }
        for(k1 = 0; k1 < NCOMP;  k1++){
          W[i][k1] = phi[i][k1] / sumi;
          sumwi[k1] += W[i][k1];
        }
      }

      /**********End of the E-step*******/

      for(i = 0; i < nsize; i++){
        for(k1 = 0; k1 < NCOMP;  k1++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * beta0[j][k1];
          mui += alpha0[k1];
          w_s[i][k1] = (resp[i] - mui) / sigma0[k1];
          Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 - pnorm(w_s[i][k1], 0, 1, 1, 0)));
          V[i][k1] = delta[i] * resp[i] + (1- delta[i]) * (mui + sigma0[k1] * Aw[i][k1]);
        }
      }

      /*********M-step of the EM*********/

      for(k1 = 0; k1 < NCOMP; k1++){

        /******Constructing the weigthed vector XTWY***/

        for(j = 0; j < NCOV + 1; j++){
          sumi1 = 0.0;
          for(i = 0; i < nsize; i++)
            sumi1 += one_X[i][j] * W[i][k1] * V[i][k1];
          oneXTWY[j] = sumi1;
        }

        /*****Constructing the weighted matrix XTWX***/

        for(i = 0 ; i < (NCOV + 1); i++){
          for(j = i; j < (NCOV + 1); j++){
            sumi2 = 0.0;
            for(l = 0; l < nsize; l++)
              sumi2 += one_X[l][i] * W[l][k1] * one_X[l][j];
            if(i == j)
              oneXTWX[i][j] = sumi2 + ridge1 * log(nsize);
            else
              oneXTWX[j][i] = oneXTWX[i][j] = sumi2;
          }
        }

        oneXTWX[0][0] = oneXTWX[0][0] - ridge1 * log(nsize);

        /***In a system Ax=b, adding b to A as its last column**/

        for(i = 0 ; i < (NCOV + 1); i++)
          for(j = 0; j < (NCOV + 2); j++)
            if(j != (NCOV + 1))
              oneComMat[i][j] = oneXTWX[i][j];
            else
              oneComMat[i][j] = oneXTWY[i];

            for (i = 0; i < (NCOV + 2); i++) {
              for (j = 0; j < (NCOV + 2); j++) {
                oneComMatVec[j + i * (NCOV + 2)] = oneComMat[i][j];
              }
            }

            /**************************************************************/
            /*Solving the system Ax=y to get betahat in the k-th component*/
            /**************************************************************/

            sol(NCOV+1, oneComMatVec, onesolution1, check1);

            for(j = 0; j < (NCOV + 1); j++)
              if(j == 0)
                new_alpha0[k1] = onesolution1[j];
              else
                new_beta0[j-1][k1] = onesolution1[j];
      }

      for(i = 0; i < nsize; i++){
        for(k1 = 0; k1 < NCOMP;  k1++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * new_beta0[j][k1];
          mui += new_alpha0[k1];
          sumi5[k1] += W[i][k1] * pow(V[i][k1] - mui, 2);
          sumi3[k1] += W[i][k1] * (delta[i] + (1 - delta[i]) * (Aw[i][k1] * (Aw[i][k1] - w_s[i][k1])))  ;
        }
      }

      for (k1 = 0; k1 < NCOMP; k1++) {
        new_sigma0[k1] = sqrt((sumi5[k1] + sigpennom) / (sumi3[k1] + sigpendenom));
      }

      for (k1 = 0; k1 < NCOMP; k1++) {
        new_pi0[k1] = sumwi[k1] / nsize;
      }

      /*****End of the M-step of the EM*****/
      jamconvg1 = 0.0;
      niter1++;
      for(k1 = 0; k1 < NCOMP;  k1++){
        for(j = 0; j < NCOV; j++)
          jamconvg1 += pow(new_beta0[j][k1] - beta0[j][k1], 2);
        jamconvg1 += pow(new_pi0[k1] - pi0[k1], 2);
        jamconvg1 += pow(new_alpha0[k1] - alpha0[k1], 2);
        jamconvg1 += pow(new_sigma0[k1] - sigma0[k1], 2);
      }

      if(jamconvg1 <= eps_conv)
        convg1='y';

      for(k1 = 0; k1 < NCOMP;  k1++){
        for(j = 0; j < NCOV; j++)
          beta0[j][k1] = new_beta0[j][k1];
        alpha0[k1] = new_alpha0[k1];
        pi0[k1] = new_pi0[k1];
        sigma0[k1] = new_sigma0[k1];
      }

      /*******End of each iteration*******/
    }

    for(k1 = 0; k1 < NCOMP;  k1++){
      for(j = 0; j < NCOV; j++)
        beta0hat[j][k1] = new_beta0[j][k1];
      alpha0hat[k1] = new_alpha0[k1];
      pi0hat[k1] = new_pi0[k1];
      sigma0hat[k1] = new_sigma0[k1];
    }

    loglike1 = 0.0;

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0hat[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0hat[k1], 1, 0), 1 - delta[i]) ;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      loglike1 += log(sumi);
    }

    *loglike = loglike1;

    for (k1 = 0; k1 < NCOMP; k1++)
    {
      for (j = 0;  j < NCOV;  j++)
      {
        mybeta[k1 * NCOV + j] = beta0hat[j][k1];
      }
      myalpha[k1] = alpha0hat[k1];
      mysigma[k1] = sigma0hat[k1];
      mypi[k1] = pi0hat[k1];
    }

    *BIC = loglike1 - 0.5 * NCOMP * NCOV * log(nsize);
    *EBIC5 = loglike1 - 0.5 * NCOMP * NCOV * log(nsize) - 0.5 * (NCOMP * NCOV) * log(NCOV);
    *EBIC1 = loglike1 - 0.5 * (NCOMP * NCOV) * log(nsize) - (NCOMP * NCOV) * log(NCOV);
    *AIC = loglike1 - (NCOMP * NCOV);
    *GCV = (loglike1) / (nsize * pow(1 - NCOMP * NCOV / nsize, 2));
    *GIC = loglike1 - 0.5 * (NCOMP * NCOV) * log(nsize);
    *MaxEMiter = niter1;

    /*******The E-step of the EM********/

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0hat[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0hat[k1], 1, 0), 1 - delta[i]) ;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      for(k1 = 0; k1 < NCOMP;  k1++){
        W[i][k1] = phi[i][k1] / sumi;
      }
    }

    /**********End of the E-step*******/
    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < nsize; i++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        predict[k1 * nsize + i] = mui;
        residual[k1 * nsize + i] = resp[i] - mui;
        tau[k1 * nsize + i] = W[i][k1];
      }
    }

  }
#ifdef __cplusplus
}
#endif

/* ********************* Variable Selection Normal and Log-Normal *********** */

#ifdef __cplusplus
extern "C" {
#endif
  void FMR_Norm_Surv_EM_VarSel(double *resp,
                               double *myX,
                               double *delta,
                               int *myPenaltyFamily,
                               double *mylambda,
                               double *myridgepen_prop,
                               int *myNCOMP,
                               int *myNCOV,
                               int *mynsize,
                               int *myEMmaxiter,
                               int *MaxEMiter,
                               double *myinitial_alpha,
                               double *myinitial_beta,
                               double *myinitial_sigma,
                               double *myinitial_pi,
                               double *myeps,
                               double *myepsconv,
                               double *gammixportion,
                               double *myalpha,
                               double *mybeta,
                               double *mysigma,
                               double *mypi,
                               double *loglike,
                               double *BIC,
                               double *AIC,
                               double *GCV,
                               double *EBIC1,
                               double *EBIC5,
                               double *GIC,
                               double *predict,
                               double *residual,
                               double *tau
  )
  {
    int nsize = *mynsize;
    int NCOV = *myNCOV;
    int EM_maxiter = *myEMmaxiter;
    int NCOMP = *myNCOMP;
    int jar = *myPenaltyFamily;

    double ridge1 = *myridgepen_prop;
    double eps = *myeps;
    double eps_conv = *myepsconv;
    double GamMP = *gammixportion;

    int i;
    int j;
    int k1;
    int l;
    int niter1;
    int check1[1];
    int SUM1;

    double multX[nsize][NCOV];
    double one_X[nsize][NCOV + 1];

    double optlam[NCOMP];

    double initbeta[NCOV][NCOMP];
    double beta0[NCOV][NCOMP];
    double new_beta0[NCOV][NCOMP];
    double beta0hat[NCOV][NCOMP];

    double initalpha[NCOMP];
    double alpha0[NCOMP];
    double new_alpha0[NCOMP];
    double alpha0hat[NCOMP];

    double initsigma[NCOMP];
    double sigma0[NCOMP];
    double new_sigma0[NCOMP];
    double sigma0hat[NCOMP];

    double sigpennom = 0.0;
    double sigpendenom = 0.0;

    double initpi[NCOMP];
    double pi0[NCOMP];
    double new_pi0[NCOMP];
    double pi0hat[NCOMP];

    double W[nsize][NCOMP];
    double phi[nsize][NCOMP];

    double w_s[nsize][NCOMP];
    double Aw[nsize][NCOMP];
    double V[nsize][NCOMP];

    double eps1[NCOMP];

    double vecder[NCOV];
    double vecsigma[NCOV + 1];
    double En[NCOMP];

    double oneXTWY[NCOV + 1];
    double oneXTWX[NCOV + 1][NCOV + 1];
    double selection[NCOV][NCOMP];

    double oneComMat[NCOV + 2][NCOV + 2];
    double onesolution1[NCOV + 2];
    double oneComMatVec[(NCOV + 2) * (NCOV + 2)];

    double sumwi[NCOMP];
    double sumi4[NCOMP];
    double sumi5[NCOMP];
    double sumi;
    double mui;
    double deni;

    char convg1 = 'n';

    double jamconvg1;
    double loglike1;
    double sat_loglike1;
    double sat_den;

    double holdveccov[NCOV];
    //    double holdveccom[NCOMP];

    for (k1 = 0; k1 < NCOMP; k1++) {
      optlam[k1]= mylambda[k1];
    }

    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < NCOV; i++){
        new_beta0[i][k1] = beta0hat[i][k1] = beta0[i][k1] = initbeta[i][k1] = myinitial_beta[NCOV * k1 + i];
      }
      new_alpha0[k1] = alpha0hat[k1] = alpha0[k1] = initalpha[k1] = myinitial_alpha[k1];
      new_sigma0[k1] = sigma0hat[k1] = sigma0[k1] = initsigma[k1] = myinitial_sigma[k1];
      new_pi0[k1] = pi0hat[k1] = pi0[k1] = initpi[k1] = myinitial_pi[k1];
    }

    for (j = 0; j < NCOV; j++){
      for (i = 0; i < nsize; i++){
        multX[i][j] = myX[nsize * j + i];
      }
    }

    for(i = 0; i < nsize; i++){
      one_X[i][0] = 1.0;
      for(j = 1; j < (NCOV + 1); j++){
        one_X[i][j] = multX[i][j - 1];
      }
    }

    for(k1 = 0; k1 < NCOMP; k1++){
      for (i = 0; i < NCOV; i++) {
        holdveccov[i] = initbeta[i][k1];
      }
      if((jar == 1) || (jar == 2))
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (2 * nsize * optlam[k1]);
      else
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (4 * nsize * optlam[k1]);
    }

    niter1=0;
    convg1 = 'n';

    while((convg1 != 'y') && (niter1 < EM_maxiter)){

      /******Beginning of each iteration******/

      /****E-step of the EM algorithm************/
      for(k1 = 0; k1 < NCOMP;  k1++){
        sumwi[k1] = 0.0;
      }

      for(i = 0; i < nsize; i++){
        sumi = 0.0;
        for(k1 = 0; k1 < NCOMP;  k1++){
          mui = 0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * beta0[j][k1];
          mui += alpha0[k1];
          deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0[k1], 1, 0), 1 - delta[i]) ;
          phi[i][k1] = pi0[k1] * deni;
          sumi += phi[i][k1];
        }
        for(k1 = 0; k1 < NCOMP;  k1++){
          W[i][k1] = phi[i][k1] / sumi;
          sumwi[k1] += W[i][k1];
        }
      }

      for (k1 = 0; k1 < NCOMP; k1++) {
        new_pi0[k1] = sumwi[k1] / nsize;
      }

      /*****End of the E-step of the EM algorithm*****/

      for(i = 0; i < nsize; i++){
        for(k1 = 0; k1 < NCOMP;  k1++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * beta0[j][k1];
          mui += alpha0[k1];
          w_s[i][k1] = (resp[i] - mui) / sigma0[k1];
          Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 - pnorm(w_s[i][k1], 0, 1, 1, 0)));
          V[i][k1] = delta[i] * resp[i] + (1- delta[i]) * (mui + sigma0[k1] * Aw[i][k1]);
        }
      }

      /****M-step of the EM algorithm*************/

      for(k1 = 0; k1 < NCOMP;  k1++){

        En[k1] = nsize * pow(pi0[k1], GamMP);

        for (j = 0; j < NCOV; j++) {
          holdveccov[j] = beta0[j][k1];
        }

        if(jar == 1)
          for(j = 0; j < NCOV; j++)
            vecder[j] = optlam[k1];
        else if (jar == 2)
          scadpen(optlam[k1], vecder, holdveccov, NCOV);
        else if (jar == 3)
          mcppen(optlam[k1], vecder, holdveccov, NCOV);
        else if (jar == 4)
          sicapen(optlam[k1], vecder, holdveccov, NCOV);
        else if (jar == 5){
          for(j = 0; j < NCOV; j++)
            vecder[j] = optlam[k1] / (fabs(initbeta[j][k1]) + eps);
        }
        else
          hardpen(optlam[k1], vecder, holdveccov, NCOV);


        for(j = 0; j < (NCOV + 1); j++)
          if(j == 0)
            vecsigma[j] = 0.0;
          else
            vecsigma[j] = vecder[j - 1] / (fabs(beta0[j - 1][k1]) + eps1[k1]);

          for(i = 0 ; i < (NCOV + 1); i++){
            for(j = 0; j < (NCOV + 1); j++){
              sumi = 0.0;
              for(l = 0; l < nsize; l++)
                sumi += one_X[l][i] * W[l][k1] * one_X[l][j];
              if(i != j)
                oneXTWX[i][j] = sumi;
              else
                oneXTWX[i][j] = sumi + En[k1] * vecsigma[j] + ridge1 * log(nsize);
            }
          }

          oneXTWX[0][0] = oneXTWX[0][0] - (En[k1] * vecsigma[0]) - ridge1 * log(nsize);

          /******Constructing the weigthed vector XTWY***/

          for(j = 0; j < (NCOV + 1); j++){
            sumi = 0.0;
            for(i = 0; i < nsize; i++)
              sumi += one_X[i][j] * W[i][k1] * V[i][k1];
            oneXTWY[j] = sumi;
          }

          /***In a system Ax=b, adding b to A as its last column**/

          for(i = 0; i < (NCOV + 1); i++)
            for(j = 0; j < (NCOV + 2); j++)
              if(j != (NCOV + 1))
                oneComMat[i][j] = oneXTWX[i][j];
              else
                oneComMat[i][j] = oneXTWY[i];

              for (i = 0; i < (NCOV + 2); i++) {
                for (j = 0; j < (NCOV + 2); j++) {
                  oneComMatVec[j + i * (NCOV + 2)] = oneComMat[i][j];
                }
              }

              /**************************************************************/
              /*Solving the system Ax=y to get betahat in the k-th component*/
              /**************************************************************/

              sol(NCOV+1, oneComMatVec, onesolution1, check1);

              for(j = 0; j < (NCOV + 1); j++){
                if(j == 0){
                  new_alpha0[k1] = onesolution1[j];
                }
                else{
                  new_beta0[j - 1][k1] = onesolution1[j];
                }
              }

      }//* End of each component

      {
        for (k1 = 0; k1 < NCOMP; k1++){
          sumi5[k1] = 0.0;
          sumi4[k1] = 0.0;
        }

        for(i = 0; i < nsize; i++){
          for(k1 = 0; k1 < NCOMP;  k1++){
            mui = 0.0;
            for(j = 0; j < NCOV; j++)
              mui += multX[i][j] * new_beta0[j][k1];
            mui += new_alpha0[k1];
            sumi5[k1] += W[i][k1] * pow(V[i][k1] - mui, 2);
          }
        }

        for (k1 = 0; k1 < NCOMP; k1++) {
          for (i = 0; i < nsize; i++) {
            sumi4[k1] += W[i][k1] * (delta[i] + (1 - delta[i]) * (Aw[i][k1] * (Aw[i][k1] - w_s[i][k1])))  ;
          }
        }

        for (k1 = 0; k1 < NCOMP; k1++) {
          new_sigma0[k1] = sqrt((sumi5[k1] + sigpennom) / (sumi4[k1] + sigpendenom));
        }
      }

      jamconvg1 = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        for(j = 0; j < NCOV; j++)
          jamconvg1 += pow(new_beta0[j][k1] - beta0[j][k1], 2);
        jamconvg1 += pow(new_pi0[k1] - pi0[k1], 2);
        jamconvg1 += pow(new_alpha0[k1] - alpha0[k1], 2);
        jamconvg1 += pow(new_sigma0[k1] - sigma0[k1], 2);
      }

      if (jamconvg1 <= eps_conv)
        convg1='y';

      niter1++;

      for(k1 = 0; k1 < NCOMP;  k1++){
        alpha0[k1] = new_alpha0[k1];
        for(j = 0; j < NCOV; j++)
          beta0[j][k1] = new_beta0[j][k1];
        pi0[k1] = new_pi0[k1];
        sigma0[k1] = new_sigma0[k1];
      }

    }//*End of each iteration

    //*********************************************************************
    //*Storing the estiamtes of the regression coefficients
    //*in a global variable called "Betahat".
    //*Selecting the finial model and storing
    //*********************************************************************

    for(k1 = 0; k1 < NCOMP;  k1++)
      for(j = 0; j < NCOV; j++)
        selection[j][k1] = 0;

    for(k1 = 0; k1 < NCOMP;  k1++){
      alpha0hat[k1] = new_alpha0[k1];
      for(j = 0; j < NCOV; j++){
        beta0hat[j][k1] = new_beta0[j][k1];
        if(fabs(beta0hat[j][k1]) <= 0.25)
          selection[j][k1] = 0;
        else
          selection[j][k1] = 1;
      }
    }
    for (k1 = 0; k1 < NCOMP; k1++) {
      pi0hat[k1] = new_pi0[k1];
    }

    for(i = 0; i < nsize; i++){
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1] * selection[j][k1];
        mui += alpha0hat[k1];
        w_s[i][k1] = (resp[i] - mui) / new_sigma0[k1];
        Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 - pnorm(w_s[i][k1], 0, 1, 1, 0)));
        V[i][k1] = delta[i] * resp[i] + (1 - delta[i]) * (mui + new_sigma0[k1] * Aw[i][k1]);
      }
    }

    for (k1 = 0; k1 < NCOMP; k1++){
      sumi5[k1] = 0.0;
      sumi4[k1] = 0.0;
    }

    for(i = 0; i < nsize; i++){
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++){
          mui +=multX[i][j] * beta0hat[j][k1] * selection[j][k1];}
        mui += alpha0hat[k1];
        sumi5[k1] += W[i][k1] * pow(V[i][k1] - mui,2);
      }
    }

    for (k1 = 0; k1 < NCOMP; k1++) {
      for (i = 0; i < nsize; i++) {
        sumi4[k1] += W[i][k1] * (delta[i] + (1 - delta[i]) * (Aw[i][k1] * (Aw[i][k1] - w_s[i][k1])))  ;
      }
    }

    for (k1 = 0; k1 < NCOMP; k1++) {
      sigma0hat[k1] = sqrt((sumi5[k1] + sigpennom) / (sumi4[k1] + sigpendenom));
    }

    loglike1 = 0.0;
    sat_loglike1 = 0.0;

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1] * selection[j][k1];
        mui += alpha0hat[k1];
        deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0hat[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0hat[k1], 1, 0), 1 - delta[i]) ;
        sat_den = pow(eps +  dnorm(0, 0, sigma0hat[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(0, 0, sigma0hat[k1], 1, 0), 1 - delta[i]) ;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      loglike1 += log(sumi);
      sat_loglike1 += log(sat_den);
    }

    SUM1 = 0;
    for(k1 = 0; k1 < NCOMP;  k1++)
      for(j = 0; j < NCOV; j++)
        SUM1 += selection[j][k1];

    *loglike = loglike1;

    for (k1 = 0; k1 < NCOMP; k1++)
    {
      for (j = 0;  j < NCOV;  j++)
      {
        mybeta[k1 * NCOV + j] = beta0hat[j][k1];
      }
      myalpha[k1] = alpha0hat[k1];
      mysigma[k1] = sigma0hat[k1];
      mypi[k1] = pi0hat[k1];
    }

    *BIC = loglike1 - 0.5 * SUM1 * log(nsize);
    *EBIC5 = loglike1 - 0.5 * (SUM1) * log(nsize) - 0.5 * (SUM1) * log(NCOV);
    *EBIC1 = loglike1 - 0.5 * (SUM1) * log(nsize) - (SUM1) * log(NCOV);
    *AIC = loglike1 - (SUM1);
    *GCV = (loglike1) / (nsize * pow(1 - SUM1 / nsize, 2));
    *GIC = loglike1 - 0.5 * (SUM1) * log(nsize);
    *MaxEMiter = niter1;
    /*******The E-step of the EM********/

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0hat[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0hat[k1], 1, 0), 1 - delta[i]) ;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      for(k1 = 0; k1 < NCOMP;  k1++){
        W[i][k1] = phi[i][k1] / sumi;
      }
    }

    /**********End of the E-step*******/
    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < nsize; i++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        predict[k1 * nsize + i] = mui;
        residual[k1 * nsize + i] = resp[i] - mui;
        tau[k1 * nsize + i] = W[i][k1];
      }
    }

  }
#ifdef __cplusplus
}
#endif

/* ****************** Tuning Parameter Normal and Log-Normal ************** */

#ifdef __cplusplus
extern "C" {
#endif
  void FMR_Norm_Surv_CwTuneParSel(double *resp,
                                  double *myX,
                                  double *delta,
                                  int *myPenaltyFamily,
                                  double *myridgepen_prop,
                                  int *myNCOMP,
                                  int *myNCOV,
                                  int *mynsize,
                                  double *myinitial_alpha,
                                  double *myinitial_beta,
                                  double *myinitial_sigma,
                                  double *myinitial_pi,
                                  double *myeps,
                                  double *myepsEM,
                                  double *gammixportion,
                                  double *optlambda
  )
  {
    int nsize = *mynsize;
    int NCOV = *myNCOV;
    int NCOMP = *myNCOMP;
    int jar = *myPenaltyFamily;

    double ridge1 = *myridgepen_prop;
    double gam = 0.0;
    if(ridge1 != 0.0)
      gam = 1;

    double eps = *myepsEM;
    double GamMP = *gammixportion;

    int i;
    int j;
    int k1;
    int l;
    int check1[1];

    int l1;
    int MaxLim = 80;

    double multX[nsize][NCOV];
    double one_X[nsize][NCOV + 1];

    double optlam[NCOMP];

    double initbeta[NCOV][NCOMP];
    double beta0[NCOV][NCOMP];
    double new_beta0[NCOV][NCOMP];

    double initalpha[NCOMP];
    double alpha0[NCOMP];
    double new_alpha0[NCOMP];

    double initsigma[NCOMP];
    double sigma0[NCOMP];
    double new_sigma0[NCOMP];

    //        double sigpennom = 0.0;
    //        double sigpendenom = 0.0;

    double initpi[NCOMP];
    double pi0[NCOMP];
    double new_pi0[NCOMP];

    double W[nsize][NCOMP];
    double phi[nsize][NCOMP];

    double w_s[nsize][NCOMP];
    double Aw[nsize][NCOMP];
    double V[nsize][NCOMP];

    double eps1[NCOMP];

    double vecder[NCOV];
    double vecsigma[NCOV + 1];
    double En[NCOMP];

    double oneXTWY[NCOV + 1];
    double oneXTWX[NCOV + 1][NCOV + 1];
    double selection[NCOV][NCOMP];

    double oneComMat[NCOV + 2][NCOV + 2];
    double onesolution1[NCOV + 2];
    double oneComMatVec[(NCOV + 2) * (NCOV + 2)];

    double sumi;
    double mui;
    double deni;

    double n1[NCOMP];
    double BIC[MaxLim][NCOMP];
    double Max_BIC[NCOMP];
    double lambda1[MaxLim];

    int count1[NCOMP][MaxLim];
    int indx1[NCOMP];

    for(l1 = 0; l1 < MaxLim; l1++)
      lambda1[l1] = 0.01 + l1 * 0.01;



    double loglike1;

    double holdveccov[NCOV];
    //    double holdveccom[NCOMP];

    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < NCOV; i++){
        new_beta0[i][k1] = beta0[i][k1] = initbeta[i][k1] = myinitial_beta[NCOV * k1 + i];
      }
      new_alpha0[k1] = alpha0[k1] = initalpha[k1] = myinitial_alpha[k1];
      new_sigma0[k1] = sigma0[k1] = initsigma[k1] = myinitial_sigma[k1];
      new_pi0[k1] = pi0[k1] = initpi[k1] = myinitial_pi[k1];
    }

    for (j = 0; j < NCOV; j++){
      for (i = 0; i < nsize; i++){
        multX[i][j] = myX[nsize * j + i];
      }
    }

    for(i = 0; i < nsize; i++){
      one_X[i][0] = 1.0;
      for(j = 1; j < (NCOV + 1); j++){
        one_X[i][j] = multX[i][j - 1];
      }
    }

    for(k1 = 0; k1 < NCOMP; k1++){
      for (i = 0; i < NCOV; i++) {
        holdveccov[i] = initbeta[i][k1];
      }
      if((jar == 1) || (jar == 2))
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (2 * nsize * optlam[k1]);
      else
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (4 * nsize * optlam[k1]);
    }


    for(k1 = 0; k1 < NCOMP; k1++)
      n1[k1] = 0.0;

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP; k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * initbeta[j][k1];
        mui += initalpha[k1];
        deni = pow(eps +  dnorm(resp[i] - mui, 0, initsigma[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, initsigma[k1], 1, 0), 1 - delta[i]) ;
        phi[i][k1] = initpi[k1] * deni;
        sumi += phi[i][k1];
      }
      for(k1 = 0; k1 < NCOMP; k1++){
        W[i][k1] = phi[i][k1] / sumi;
        //n1[k1] += W[i][k1];
      }
    }

    //*Start of choosing lambda for each component of the mixture

    for(k1 = 0; k1 < NCOMP; k1++){
      En[k1] = nsize * pow(initpi[k1], GamMP);

      for (j = 0; j < NCOV; j++) {
        holdveccov[j] = beta0[j][k1];
      }

      for(l1 = 0; l1 < MaxLim; l1++){
        if(jar == 1)
          for(j = 0; j < NCOV; j++)
            vecder[j] = lambda1[l1];
        else if (jar == 2)
          scadpen(lambda1[l1], vecder, holdveccov, NCOV);
        else if (jar == 3)
          mcppen(lambda1[l1], vecder, holdveccov, NCOV);
        else if (jar == 4)
          sicapen(lambda1[l1], vecder, holdveccov, NCOV);
        else if (jar == 5){
          for(j = 0; j < NCOV; j++)
            vecder[j] = optlam[k1] / (fabs(initbeta[j][k1]) + eps);
        }
        else
          hardpen(optlam[k1], vecder, holdveccov, NCOV);


        for(j = 0; j <(NCOV + 1); j++)
          if(j == 0)
            vecsigma[j] = 0.0;
          else
            vecsigma[j] = vecder[j - 1] / (fabs(beta0[j - 1][k1]) + eps1[k1]);

          for(i = 0; i < (NCOV + 1); i++){
            for(j = 0; j < (NCOV + 1); j++){
              sumi = 0.0;
              for(l = 0; l < nsize; l++)
                sumi += one_X[l][i] * W[l][k1] * one_X[l][j];
              if(i != j)
                oneXTWX[i][j] = sumi;
              else
                oneXTWX[i][j] = sumi + En[k1] * vecsigma[j] + ridge1 * log(nsize);
            }
          }

          oneXTWX[0][0] = oneXTWX[0][0] - (En[k1] * vecsigma[0]) - ridge1 * log(nsize);

          for(i = 0; i < nsize; i++){
            mui = 0.0;
            for(j = 0; j < NCOV; j++)
              mui += multX[i][j] * beta0[j][k1];
            mui += alpha0[k1];
            w_s[i][k1] = (resp[i] - mui) / sigma0[k1];
            Aw[i][k1] = dnorm(w_s[i][k1], 0, 1, 0) / (eps + (1 - pnorm(w_s[i][k1], 0, 1, 1, 0)));
            V[i][k1] = delta[i] * resp[i] + (1 - delta[i]) * (mui + sigma0[k1] * Aw[i][k1]);
          }
          /****** Constructing the weigthed vector oneXTWY ***/

          for(j = 0; j < (NCOV + 1); j++){
            sumi = 0.0;
            for(i = 0; i < nsize; i++)
              sumi += one_X[i][j] * W[i][k1] * V[i][k1];
            oneXTWY[j] = sumi;
          }
          /***In a system Ax=b, adding b to A as its last column**/

          for(i = 0; i < (NCOV + 1); i++)
            for(j = 0; j < (NCOV + 2); j++)
              if(j != (NCOV + 1))
                oneComMat[i][j] = oneXTWX[i][j];
              else
                oneComMat[i][j] = oneXTWY[i];

              for (i = 0; i < (NCOV + 2); i++) {
                for (j = 0; j < (NCOV + 2); j++) {
                  oneComMatVec[j + i * (NCOV + 2)] = oneComMat[i][j];
                }
              }

              /**************************************************************/
              /*Solving the system Ax=y to get betahat in the k-th component*/
              /**************************************************************/

              count1[k1][l1] = 0;
              sol(NCOV + 1, oneComMatVec, onesolution1, check1);

              for(j = 0; j < (NCOV + 1); j++){
                if(j == 0)
                  new_alpha0[k1] = onesolution1[j];
                else {
                  new_beta0[j - 1][k1]  = onesolution1[j];
                  if(fabs(new_beta0[j - 1][k1]) < 0.2)
                    selection[j - 1][k1] = 0;
                  else
                    selection[j - 1][k1] = 1;
                  count1[k1][l1] += selection[j - 1][k1];
                }
              }

              for(j = 0; j < NCOV; j++){
                beta0[j][k1] = new_beta0[j][k1];
              }

              alpha0[k1] = new_alpha0[k1];

              loglike1 = 0.0;
              n1[k1] = 0.0;

              for(i = 0; i < nsize; i++){
                mui = 0.0;
                for(j = 0; j < NCOV; j++)
                  mui += multX[i][j] * new_beta0[j][k1] * selection[j][k1];
                mui += alpha0[k1];
                deni = pow(eps +  dnorm(resp[i] - mui, 0, sigma0[k1], 0), delta[i]) *  pow(eps + 1 - pnorm(resp[i] - mui, 0, sigma0[k1], 1, 0), 1 - delta[i]) ;
                loglike1 += W[i][k1] * log(deni);
                n1[k1] += W[i][k1];
              }

              BIC[l1][k1] = loglike1 - 0.5 * (count1[k1][l1]) * log(n1[k1]) - gam * (count1[k1][l1]) * log(NCOV);

              if(l1 == 0){
                Max_BIC[k1] = BIC[l1][k1];
                indx1[k1] = l1;
              }
              else if(BIC[l1][k1] > Max_BIC[k1])
              {
                Max_BIC[k1] = BIC[l1][k1];
                indx1[k1] = l1;
              }
      }//*End of choosing lambda for each component
    }//*End of choosing lambda for both components of the mixture

    for(k1 = 0; k1 < NCOMP; k1++){
      optlambda[k1] = lambda1[indx1[k1]];
    }
  }
#ifdef __cplusplus
}
#endif

/* *************************** MLE Weibull *************************** */

#ifdef __cplusplus
extern "C" {
#endif
  void FMR_Weibl_Surv_EM_MLE(double *resp,
                             double *myX,
                             double *delta,
                             double *myridgepen_prop,
                             int *myNCOMP,
                             int *myNCOV,
                             int *mynsize,
                             int *myEMmaxiter,
                             int *myNRmaxiter,
                             int *myNRportion,
                             int *MaxEMiter,
                             double *myinitial_alpha,
                             double *myinitial_beta,
                             double *myinitial_sigma,
                             double *myinitial_pi,
                             double *myepsconv,
                             double *myalpha,
                             double *mybeta,
                             double *mysigma,
                             double *mypi,
                             double *loglike,
                             double *BIC,
                             double *AIC,
                             double *GCV,
                             double *EBIC1,
                             double *EBIC5,
                             double *GIC,
                             double *predict,
                             double *residual,
                             double *tau
  )
  {
    int nsize = *mynsize;
    int NCOV = *myNCOV;
    int NCOMP = *myNCOMP;
    int EM_maxiter = *myEMmaxiter;
    int NR_maxiter = *myNRmaxiter;

    double eps_conv = *myepsconv;

    double ridge1 = *myridgepen_prop;
    int NRportion = *myNRportion;
    int alp = 0;

    int i;
    int j;
    int j1;
    int k1;
    int l;
    int niter1;
    int niter2;
    int check1[1];

    double multX[nsize][NCOV];
    double one_X[nsize][NCOV + 1];
    double sumi;
    double mui;
    double deni;

    double W[nsize][NCOMP];
    double phi[nsize][NCOMP];

    double initbeta[NCOV][NCOMP];
    double beta0[NCOV][NCOMP];
    double new_beta0[NCOV][NCOMP];
    double beta0hat[NCOV][NCOMP];

    double initalpha[NCOMP];
    double alpha0[NCOMP];
    double new_alpha0[NCOMP];
    double alpha0hat[NCOMP];

    double initsigma[NCOMP];
    double sigma0[NCOMP];
    double new_sigma0[NCOMP];
    double sigma0hat[NCOMP];

    //        double sigpennom = 0.0;
    //        double sigpendenom = 0.0;

    double initpi[NCOMP];
    double pi0[NCOMP];
    double new_pi0[NCOMP];
    double pi0hat[NCOMP];

    double oneXTWY[NCOV + 1];
    double oneXTWX[NCOV + 1][NCOV + 1];

    double sumwi[NCOMP];
    double sumi3[NCOMP];
    double sumi5[NCOMP];

    double jamconvg1;

    double oneComMat[NCOV + 2][NCOV + 2];
    double onesolution1[NCOV + 2];
    double oneComMatVec[(NCOV + 2) * (NCOV + 2)];
    double loglike1;
    double oldloglike1;
    double newloglike1;

    char convg1 = 'n';
    char convg2 = 'n';

    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < NCOV; i++){
        new_beta0[i][k1] = beta0hat[i][k1] = beta0[i][k1] = initbeta[i][k1] = myinitial_beta[NCOV * k1 + i];
      }
      new_alpha0[k1] = alpha0hat[k1] = alpha0[k1] = initalpha[k1] = myinitial_alpha[k1];
      new_sigma0[k1] = sigma0hat[k1] = sigma0[k1] = initsigma[k1] = myinitial_sigma[k1];
      new_pi0[k1] = pi0hat[k1] = pi0[k1] = initpi[k1] = myinitial_pi[k1];
    }

    for (j = 0; j < NCOV; j++){
      for (i = 0; i < nsize; i++){
        multX[i][j] = myX[nsize * j + i];
      }
    }

    for(i = 0; i < nsize; i++){
      one_X[i][0] = 1.0;
      for(j = 1; j < (NCOV + 1); j++){
        one_X[i][j] = multX[i][j - 1];
      }
    }

    niter1 = 0;
    convg1 = 'n';
    while((convg1 != 'y') && (niter1 < EM_maxiter)){

      /****Beginning of each iteration****/

      /*******The E-step of the EM********/
      for (k1 = 0; k1 < NCOMP; k1++) {
        sumwi[k1] = 0.0;
      }

      for(i = 0; i < nsize; i++){
        sumi = 0.0;
        for(k1 = 0; k1 < NCOMP; k1++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * beta0[j][k1];
          mui += alpha0[k1];
          deni = pow((1 / sigma0[k1])* exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(-exp((resp[i] - mui) / sigma0[k1])) ;
          //if (isnan(deni))
          //  cout << "NAN" << "\t";
          //                    if (deni < eps2)
          //                        deni = eps2;
          //cout << deni << endl;
          phi[i][k1] = pi0[k1] * deni;
          sumi += phi[i][k1];
        }
        for(k1 = 0; k1 < NCOMP; k1++){
          W[i][k1] = phi[i][k1] / sumi;
          sumwi[k1] += W[i][k1];
        }
      }

      for(k1 = 0; k1 < NCOMP; k1++){
        new_pi0[k1] = sumwi[k1] / nsize;
      }

      /**********End of the E-step*******/

      for (k1 = 0; k1 < NCOMP; k1++) {
        for (j1 = 0; j1 < NCOV; j1++) {
          initbeta[j1][k1] = beta0[j1][k1];
        }
        initalpha[k1] = alpha0[k1];
        initsigma[k1] = sigma0[k1];
        initpi[k1] = pi0[k1];
      }

      /*********M-step of the EM*********/

      for(k1 = 0; k1 < NCOMP; k1++){
        niter2 = 0;
        convg2 = 'n';
        while((convg2 != 'y') && (niter2 < NR_maxiter)){

          /******Constructing the weigthed score function***/
          oldloglike1 = 0.0;
          for(i = 0; i < nsize; i++){
            mui = 0.0;
            for(j = 0; j < NCOV; j++)
              mui += multX[i][j] * beta0[j][k1];
            mui += alpha0[k1];
            deni = pow((1 / sigma0[k1]) * exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0[k1])) ;
            //                        if (deni < eps2)
            //                            deni = eps2;
            oldloglike1 += W[i][k1] * log(deni);
          }

          for(j = 0; j < (NCOV + 1); j++){
            sumi = 0.0;
            for(i = 0; i < nsize; i++){
              mui = 0.0;
              for(l = 0; l < NCOV; l++)
                mui += multX[i][l] * beta0[l][k1];
              mui += alpha0[k1];
              sumi += W[i][k1] * one_X[i][j] / sigma0[k1] * (exp((resp[i] - mui) / sigma0[k1]) - delta[i]);
            }
            oneXTWY[j] =  sumi;
          }

          /*****Constructing the weighted hessian matrix***/

          for(i = 0; i < (NCOV + 1); i++){
            for(j = i; j < (NCOV + 1); j++){
              sumi = 0.0;
              for(l = 0; l < nsize; l++){
                mui = 0.0;
                for(j1 = 0; j1 < NCOV; j1++)
                  mui += multX[l][j1] * beta0[j1][k1];
                mui += alpha0[k1];
                sumi += - W[l][k1] * one_X[l][i] * one_X[l][j] / (sigma0[k1] * sigma0hat[k1]) * exp((resp[l] - mui) / sigma0[k1]);
              }
              if(i == j)
                oneXTWX[i][j] = sumi + ridge1 * log(nsize);
              else
                oneXTWX[j][i] = oneXTWX[i][j] = sumi;
            }
          }

          oneXTWX[0][0] = oneXTWX[0][0] - ridge1 * log(nsize);

          /***In a system Ax=b, adding b to A as its last column**/

          for(i = 0; i < (NCOV + 1); i++)
            for(j = 0; j < (NCOV + 2); j++)
              if(j != (NCOV + 1))
                oneComMat[i][j] = - oneXTWX[i][j];
              else
                oneComMat[i][j] = oneXTWY[i];

              for (i = 0; i < (NCOV + 2); i++) {
                for (j = 0; j < (NCOV + 2); j++) {
                  oneComMatVec[j + i * (NCOV + 2)] = oneComMat[i][j];
                }
              }

              /**************************************************************/
              /*Solving the system Ax=y to get betahat in the k-th component*/
              /**************************************************************/

              sol((NCOV + 1), oneComMatVec, onesolution1, check1);

              alp = 0;
              do {
                for(j = 0; j < (NCOV + 1); j++)
                  if(j == 0)
                    new_alpha0[k1] = pow(0.5, alp) * onesolution1[j] + alpha0[k1];
                  else
                    new_beta0[j-1][k1] = pow(0.5, alp) * onesolution1[j] + beta0[j-1][k1];

                  sumi3[k1] = 0.0;
                  sumi5[k1] = 0.0;

                  for(i = 0; i < nsize; i++){
                    mui = 0.0;
                    for(l = 0; l < NCOV; l++)
                      mui += multX[i][l] * new_beta0[l][k1];
                    mui += new_alpha0[k1];
                    sumi3[k1] += W[i][k1] * (- delta[i] / sigma0[k1] + ((resp[i] - mui) / (sigma0[k1] * sigma0[k1])) * ( exp( (resp[i] - mui) / sigma0[k1]) - delta[i])   );
                    sumi5[k1] += W[i][k1] * ( delta[i] / (sigma0[k1] * sigma0[k1]) +  ((resp[i] - mui) / (sigma0[k1] * sigma0[k1] * sigma0[k1])) * (2 * delta[i] - (2 + (resp[i] - mui) / sigma0[k1]) * exp( (resp[i] - mui) / sigma0[k1])  )) ;
                  }
                  //sumi5[k1] += sigpennom;
                  //sumi3[k1] += sigpendenom;
                  new_sigma0[k1] = sigma0[k1] - pow(0.5, alp) * (1 / sumi5[k1]) * sumi3[k1];
                  //                        new_sigma0[k1] = (new_sigma0[k1] < 0.1)?0.5:new_sigma0[k1];
                  //                        if (k1 == 0)
                  //                            new_sigma0[k1] = (new_sigma0[k1] > 10)?2:(new_sigma0[k1]); // + 0.01
                  //                        else
                  //                            newsigma[k1] = (new_sigma0[k1] > 10)?2:(new_sigma0[k1]); // + 0.005

                  newloglike1 = 0.0;
                  for(i = 0; i < nsize; i++){
                    mui = 0.0;
                    for(j = 0; j < NCOV; j++)
                      mui += multX[i][j] * new_beta0[j][k1];
                    mui += new_alpha0[k1];
                    deni = pow((1 / new_sigma0[k1]) * exp((resp[i] - mui) / new_sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / new_sigma0[k1])) ;
                    //                            if (deni < eps2)
                    //                                deni = eps2;
                    newloglike1 += W[i][k1] * log(deni);
                  }
                  alp++;

              } while ((oldloglike1 > newloglike1) & (alp < NRportion));

              jamconvg1 = 0.0;
              niter2++;
              for(j = 0; j < NCOV; j++)
                jamconvg1 += pow(new_beta0[j][k1] - beta0[j][k1], 2);
              jamconvg1 += pow(new_alpha0[k1] - alpha0[k1], 2);
              jamconvg1 += pow(new_sigma0[k1] - sigma0[k1], 2);
              jamconvg1 += pow(new_pi0[k1] - pi0[k1], 2);

              convg2 = 'n';
              if(jamconvg1 < eps_conv)
                convg2 = 'y';

              for(j = 0; j < NCOV; j++)
                beta0[j][k1] = new_beta0[j][k1];
              alpha0[k1] = new_alpha0[k1];
              sigma0[k1] = new_sigma0[k1];
        }
      }

      /*****End of the M-step of the EM*****/

      jamconvg1 = 0.0;
      niter1++;
      for(k1 = 0; k1 < NCOMP; k1++){
        for(j = 0; j < NCOV; j++)
          jamconvg1 += pow(new_beta0[j][k1] - initbeta[j][k1], 2);
        jamconvg1 += pow(new_alpha0[k1] - initalpha[k1], 2);
        jamconvg1 += pow(new_sigma0[k1] - initsigma[k1], 2);
        jamconvg1 += pow(new_pi0[k1] - initpi[k1], 2);
      }

      convg1 = 'n';
      if(jamconvg1 < eps_conv)
        convg1 = 'y';
      for(k1 = 0; k1 < NCOMP; k1++){
        alpha0[k1] = new_alpha0[k1];
        for(j = 0; j < NCOV; j++){
          beta0[j][k1] = new_beta0[j][k1];
        }
        sigma0[k1] = new_sigma0[k1];
        pi0[k1] = new_pi0[k1];
      }
    }

    /*******End of each iteration *******/

    for(k1 = 0; k1 < NCOMP; k1++){
      for(j = 0; j < NCOV; j++){
        beta0hat[j][k1] = new_beta0[j][k1];
      }
      alpha0hat[k1] = new_alpha0[k1];
      pi0hat[k1] = new_pi0[k1];
      sigma0hat[k1] = new_sigma0[k1];
    }

    loglike1 = 0.0;

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP; k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        deni = pow((1 / sigma0hat[k1]) * exp((resp[i] - mui) / sigma0hat[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0hat[k1]));
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      loglike1 += log(sumi);
    }

    *loglike = loglike1;

    for (k1 = 0; k1 < NCOMP; k1++)
    {
      for (j = 0;  j < NCOV;  j++)
      {
        mybeta[k1 * NCOV + j] = beta0hat[j][k1];
      }
      myalpha[k1] = alpha0hat[k1];
      mysigma[k1] = sigma0hat[k1];
      mypi[k1] = pi0hat[k1];
    }

    *BIC = loglike1 - 0.5 * NCOMP * NCOV * log(nsize);
    *EBIC5 = loglike1 - 0.5 * NCOMP * NCOV * log(nsize) - 0.5 * (NCOMP * NCOV) * log(NCOV);
    *EBIC1 = loglike1 - 0.5 * (NCOMP * NCOV) * log(nsize) - (NCOMP * NCOV) * log(NCOV);
    *AIC = loglike1 - (NCOMP * NCOV);
    *GCV = (loglike1) / (nsize * pow(1 - NCOMP * NCOV / nsize, 2));
    *GIC = loglike1 - 0.5 * (NCOMP * NCOV) * log(nsize);
    *MaxEMiter = niter1;

    /*******The E-step of the EM********/

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        deni = pow((1 / sigma0hat[k1])* exp((resp[i] - mui) / sigma0hat[k1]), delta[i]) * exp(-exp((resp[i] - mui) / sigma0hat[k1])) ;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      for(k1 = 0; k1 < NCOMP;  k1++){
        W[i][k1] = phi[i][k1] / sumi;
      }
    }

    /**********End of the E-step*******/
    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < nsize; i++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        predict[k1 * nsize + i] = mui;
        residual[k1 * nsize + i] = resp[i] - mui;
        tau[k1 * nsize + i] = W[i][k1];
      }
    }

  }
#ifdef __cplusplus
}
#endif

/* ************************* Variable Selection Weibull***************** */

#ifdef __cplusplus
extern "C" {
#endif
  void FMR_Weibl_Surv_EM_VarSel(double *resp,
                                double *myX,
                                double *delta,
                                int *myPenaltyFamily,
                                double *mylambda,
                                double *myridgepen_prop,
                                int *myNCOMP,
                                int *myNCOV,
                                int *mynsize,
                                int *myEMmaxiter,
                                int *myNRmaxiter,
                                int *myNRportion,
                                int *MaxEMiter,
                                double *myinitial_alpha,
                                double *myinitial_beta,
                                double *myinitial_sigma,
                                double *myinitial_pi,
                                double *myeps,
                                double *myepsconv,
                                double *gammixportion,
                                double *myalpha,
                                double *mybeta,
                                double *mysigma,
                                double *mypi,
                                double *loglike,
                                double *BIC,
                                double *AIC,
                                double *GCV,
                                double *EBIC1,
                                double *EBIC5,
                                double *GIC,
                                double *predict,
                                double *residual,
                                double *tau
  )
  {
    int nsize = *mynsize;
    int NCOV = *myNCOV;
    int NCOMP = *myNCOMP;
    int EM_maxiter = *myEMmaxiter;
    int NR_maxiter = *myNRmaxiter;
    int NRportion = *myNRportion;
    int jar = *myPenaltyFamily;

    double ridge1 = *myridgepen_prop;
    double eps = *myeps;
    double eps_conv = *myepsconv;
    double GamMP = *gammixportion;

    int i;
    int j;
    int j1;
    int k1;
    int l;
    int niter1;
    int niter2;
    int check1[1];
    int SUM1;

    int alp = 0;

    double multX[nsize][NCOV];
    double one_X[nsize][NCOV + 1];

    double optlam[NCOMP];

    double initbeta[NCOV][NCOMP];
    double beta0[NCOV][NCOMP];
    double new_beta0[NCOV][NCOMP];
    double beta0hat[NCOV][NCOMP];

    double initalpha[NCOMP];
    double alpha0[NCOMP];
    double new_alpha0[NCOMP];
    double alpha0hat[NCOMP];

    double initsigma[NCOMP];
    double sigma0[NCOMP];
    double new_sigma0[NCOMP];
    double sigma0hat[NCOMP];

    double initpi[NCOMP];
    double pi0[NCOMP];
    double new_pi0[NCOMP];
    double pi0hat[NCOMP];

    double W[nsize][NCOMP];
    double phi[nsize][NCOMP];

    double eps1[NCOMP];

    double vecder[NCOV];
    double vecsigma[NCOV + 1];
    double En[NCOMP];

    double oneXTWY[NCOV + 1];
    double oneXTWX[NCOV + 1][NCOV + 1];
    double selection[NCOV][NCOMP];

    double oneComMat[NCOV + 2][NCOV + 2];
    double onesolution1[NCOV + 2];
    double oneComMatVec[(NCOV + 2) * (NCOV + 2)];

    double sumwi[NCOMP];
    double sumi3[NCOMP];
    double sumi5[NCOMP];
    double sumi;
    double mui;
    double deni;

    char convg1 = 'n';
    char convg2 = 'n';

    double jamconvg1;
    double loglike1;
    double oldloglike1;
    double newloglike1;
    double sat_loglike1;
    double sat_den;

    double holdveccov[NCOV];
    //    double holdveccom[NCOMP];

    for (k1 = 0; k1 < NCOMP; k1++) {
      optlam[k1]= mylambda[k1];
    }

    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < NCOV; i++){
        new_beta0[i][k1] = beta0hat[i][k1] = beta0[i][k1] = initbeta[i][k1] = myinitial_beta[NCOV * k1 + i];
      }
      new_alpha0[k1] = alpha0hat[k1] = alpha0[k1] = initalpha[k1] = myinitial_alpha[k1];
      new_sigma0[k1] = sigma0hat[k1] = sigma0[k1] = initsigma[k1] = myinitial_sigma[k1];
      new_pi0[k1] = pi0hat[k1] = pi0[k1] = initpi[k1] = myinitial_pi[k1];
    }

    for (j = 0; j < NCOV; j++){
      for (i = 0; i < nsize; i++){
        multX[i][j] = myX[nsize * j + i];
      }
    }

    for(i = 0; i < nsize; i++){
      one_X[i][0] = 1.0;
      for(j = 1; j < (NCOV + 1); j++){
        one_X[i][j] = multX[i][j - 1];
      }
    }

    for(k1 = 0; k1 < NCOMP; k1++){
      for (i = 0; i < NCOV; i++) {
        holdveccov[i] = initbeta[i][k1];
      }
      if((jar == 1) || (jar == 2))
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (2 * nsize * optlam[k1]);
      else
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (4 * nsize * optlam[k1]);
    }

    niter1 = 0;
    convg1 = 'n';
    while((convg1 != 'y') && (niter1 < EM_maxiter)){

      /******Beginning of each iteration******/

      /****E-step of the EM algorithm************/

      for(k1 = 0; k1 < NCOMP; k1++){
        sumwi[k1] = 0.0;
      }

      for(i = 0; i < nsize; i++){
        sumi = 0.0;
        for(k1 = 0; k1 < NCOMP; k1++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * beta0[j][k1];
          mui += alpha0[k1];
          deni = pow((1 / sigma0[k1]) * exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0[k1])) ;
          //                    if (isnan(deni))
          //                        cout << "NAN" << "\t";
          //                    if (deni < eps2)
          //                        deni = eps2;
          phi[i][k1] = pi0[k1] * deni;
          sumi += phi[i][k1];

        }
        for(k1 = 0; k1 < NCOMP; k1++){
          W[i][k1] = phi[i][k1] / sumi;
          sumwi[k1] += W[i][k1];
        }
      }

      for(k1 = 0; k1 < NCOMP; k1++){
        new_pi0[k1] = sumwi[k1] / nsize;
      }

      /*****End of the E-step of the EM algorithm*****/

      for (k1 = 0; k1 < NCOMP; k1++) {
        for (j1 = 0; j1 < NCOV; j1++) {
          initbeta[j1][k1] = beta0[j1][k1];
        }
        initalpha[k1] = alpha0[k1];
        initsigma[k1] = sigma0[k1];
        initpi[k1] = pi0[k1];
      }

      /****M-step of the EM algorithm*************/

      for(k1 = 0; k1 < NCOMP; k1++){
        niter2 = 0;
        convg2 = 'n';
        while((convg2 != 'y') && (niter2 < NR_maxiter)){

          oldloglike1 = 0.0;
          for(i = 0; i < nsize; i++){
            mui = 0.0;
            for(j = 0; j < NCOV; j++)
              mui += multX[i][j] * beta0[j][k1];
            mui += alpha0[k1];
            deni = pow((1 / sigma0[k1]) * exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0[k1])) ;
            //                    if (deni<eps2)
            //                        deni=eps2;
            oldloglike1 += W[i][k1] * log(deni);
          }

          En[k1] = nsize * pow(pi0[k1], GamMP);


          for (j = 0; j < NCOV; j++) {
            holdveccov[j] = beta0[j][k1];
          }

          if(jar == 1)
            for(j = 0; j < NCOV; j++)
              vecder[j] = optlam[k1];
          else if (jar == 2)
            scadpen(optlam[k1], vecder, holdveccov, NCOV);
          else if (jar == 3)
            mcppen(optlam[k1], vecder, holdveccov, NCOV);
          else if (jar == 4)
            sicapen(optlam[k1], vecder, holdveccov, NCOV);
          else if (jar == 5){
            for(j = 0; j < NCOV; j++)
              vecder[j] = optlam[k1] / (fabs(initbeta[j][k1]) + eps);
          }
          else
            hardpen(optlam[k1], vecder, holdveccov, NCOV);


          for(j = 0; j < (NCOV + 1); j++)
            if(j == 0)
              vecsigma[j] = 0.0;
            else
              vecsigma[j] = vecder[j-1] / (fabs(beta0[j - 1][k1]) + eps1[k1]);

            /******Constructing the Hessian matrix H ***/
            for(i = 0; i < (NCOV + 1); i++){
              for(j = i; j < (NCOV + 1); j++){
                sumi = 0.0;
                for(l = 0; l < nsize; l++){
                  mui = 0.0;
                  for(j1 = 0; j1 < NCOV; j1++)
                    mui += multX[l][j1] * beta0[j1][k1];
                  mui += alpha0[k1];
                  sumi += - W[l][k1] * one_X[l][i] * one_X[l][j] / (sigma0[k1] * sigma0[k1]) * exp((resp[l] - mui) / sigma0[k1]);
                }
                if(i == j)
                  oneXTWX[i][j] = sumi - En[k1] * vecsigma[j] + ridge1 * log(nsize);
                else
                  oneXTWX[j][i] = oneXTWX[i][j] = sumi;
              }
            }

            oneXTWX[0][0] = oneXTWX[0][0] + (En[k1] * vecsigma[0]) - ridge1 * log(nsize);

            /******Constructing the weigthed vector XTWY***/

            for(j = 0; j < (NCOV + 1); j++){
              sumi = 0.0;
              for(i = 0; i < nsize; i++){
                mui = 0.0;
                for(l = 0; l < NCOV; l++)
                  mui += multX[i][l] * beta0[l][k1];
                mui += alpha0[k1];
                sumi += W[i][k1] * (one_X[i][j] / sigma0[k1]) * (exp((resp[i] - mui) / sigma0[k1]) - delta[i]);
              }
              if (j == 0)
                oneXTWY[j] =  sumi - En[k1] * vecsigma[j];
              else
                oneXTWY[j] =  sumi - En[k1] * vecsigma[j] * beta0[j-1][k1];
            }

            /***In a system Ax=b, adding b to A as its last column**/

            for(i = 0; i < (NCOV + 1); i++)
              for(j = 0; j < (NCOV + 2); j++)
                if(j != (NCOV + 1))
                  oneComMat[i][j] = - oneXTWX[i][j];
                else
                  oneComMat[i][j] = oneXTWY[i];

                for (i = 0; i < (NCOV + 2); i++) {
                  for (j = 0; j < (NCOV + 2); j++) {
                    oneComMatVec[j + i * (NCOV + 2)] = oneComMat[i][j];
                  }
                }

                /**************************************************************/
                /*Solving the system Ax=y to get betahat in the k-th component*/
                /**************************************************************/

                sol( (NCOV + 1) , oneComMatVec, onesolution1, check1);

                alp = 0;
                do {
                  for(j = 0; j < (NCOV + 1); j++)
                    if(j == 0)
                      new_alpha0[k1] = pow(0.5, alp) * onesolution1[j] + alpha0[k1];
                    else
                      new_beta0[j-1][k1] = pow(0.5, alp) * onesolution1[j] + beta0[j-1][k1];

                    sumi5[k1] = 0.0;
                    sumi3[k1] = 0.0;

                    for(i = 0; i < nsize; i++){
                      mui = 0.0;
                      for(l = 0; l < NCOV; l++)
                        mui += multX[i][l] * new_beta0[l][k1];
                      mui += new_alpha0[k1];
                      sumi3[k1] += W[i][k1] * (- delta[i] / sigma0[k1] + ((resp[i] - mui) / (sigma0[k1] * sigma0[k1])) * ( exp( (resp[i] - mui) / sigma0[k1]) - delta[i])   );
                      sumi5[k1] += W[i][k1] * ( delta[i] / (sigma0[k1] * sigma0[k1]) +  ((resp[i] - mui) / (sigma0[k1] * sigma0[k1] * sigma0[k1])) * (2 * delta[i] - (2 + (resp[i] - mui) / sigma0[k1]) * exp( (resp[i] - mui) / sigma0[k1])  )) ;
                    }
                    //sumi5[k1]+=sigpennom;
                    //sumi3[k1]+=sigpendenom;
                    new_sigma0[k1] = sigma0[k1] - pow(0.5, alp) * (1 / sumi5[k1]) * sumi3[k1];
                    //                    new_sigma0[k1] = (new_sigma0[k1] < 0.6)?0.6:new_sigma0[k1];
                    //                    new_sigma0[k1] = (new_sigma0[k1] > 4)?2:new_sigma0[k1];
                    //                    if (k1==0)
                    //                        new_sigma0[k1] =  (new_sigma0[k1]>4)?2:(new_sigma0[k1]);
                    //                    else
                    //                        new_sigma0[k1] =  (new_sigma0[k1]>4)?2:(new_sigma0[k1]);

                    newloglike1 = 0.0;
                    for(i = 0; i < nsize; i++){
                      mui = 0.0;
                      for(j = 0; j < NCOV; j++)
                        mui += multX[i][j] * new_beta0[j][k1];
                      mui += new_alpha0[k1];
                      deni = pow((1 / new_sigma0[k1]) * exp((resp[i] - mui) / new_sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / new_sigma0[k1])) ;
                      //                        if (deni < eps2)
                      //                            deni = eps2;
                      newloglike1 += W[i][k1] * log(deni);
                    }
                    alp++;

                } while ((oldloglike1 > newloglike1) & (alp < NRportion));

                jamconvg1 = 0.0;
                niter2++;

                for(j = 0; j < NCOV; j++)
                  jamconvg1 += pow(new_beta0[j][k1] - beta0[j][k1], 2);
                jamconvg1 += pow(new_pi0[k1] - pi0[k1], 2);
                jamconvg1 += pow(new_alpha0[k1] - alpha0[k1], 2);
                jamconvg1 += pow(new_sigma0[k1] - sigma0[k1], 2);

                convg2 = 'n';
                if (jamconvg1 < eps)
                  convg2 = 'y';

                for(j = 0; j < NCOV; j++)
                  beta0[j][k1] = new_beta0[j][k1];
                alpha0[k1] = new_alpha0[k1];
                sigma0[k1] = new_sigma0[k1];
        }

      }
      /*****End of the M-step of the EM*****/

      jamconvg1 = 0.0;
      for(k1 = 0; k1 < NCOMP; k1++){
        for(j = 0; j < NCOV; j++)
          jamconvg1 += pow(new_beta0[j][k1] - initbeta[j][k1], 2);
        jamconvg1 += pow(new_alpha0[k1] - initalpha[k1], 2);
        jamconvg1 += pow(new_sigma0[k1] - initsigma[k1], 2);
        jamconvg1 += pow(new_pi0[k1] - initpi[k1], 2);
      }

      convg1 = 'n';
      if(jamconvg1 < eps_conv)
        convg1 = 'y';
      niter1++;

      for(k1 = 0; k1 < NCOMP; k1++){
        alpha0[k1] = new_alpha0[k1];
        for(j = 0; j < NCOV; j++){
          beta0[j][k1] = new_beta0[j][k1];
        }
        sigma0[k1] = new_sigma0[k1];
        pi0[k1] = new_pi0[k1];
      }
    }
    /*******End of each iteration *******/


    //*********************************************************************
    //*Storing the estiamtes of the regression coefficients
    //*in a global variable called "Betahat".
    //*Selecting the finial model and storing
    //*********************************************************************

    for(k1 = 0; k1 < NCOMP; k1++)
      for(j = 0; j < NCOV; j++)
        selection[j][k1] = 0;

    for(k1 = 0; k1 < NCOMP; k1++){
      alpha0hat[k1] = new_alpha0[k1];
      sigma0hat[k1] = sigma0[k1];
      pi0hat[k1] = pi0[k1];
      pi0hat[k1] = new_pi0[k1];
      for(j = 0; j < NCOV; j++){
        beta0hat[j][k1] = new_beta0[j][k1];
        if(fabs(beta0hat[j][k1]) <= 0.1)
          selection[j][k1] = 0;
        else
          selection[j][k1] = 1;
      }
    }

    //        for (k1 = 0; k1 < NCOMP; k1++){
    //            sumi5[k1] = 0.0;
    //            sumi3[k1] = 0.0;
    //
    //            for(i = 0; i < nsize; i++){
    //                mui = 0.0;
    //                for(l = 0; l < NCOV; l++)
    //                    mui += multX[i][l] * beta0hat[l][k1] * selection[j][k1];
    //                mui += alpha0hat[k1];
    //                sumi3[k1] += W[i][k1] * (- 1 / sigma0[k1] + ((resp[i] - mui) / (sigma0[k1] * new_sigma0[k1])) * ( exp( (resp[i] - mui) / sigma0[k1]) - delta[i])   );
    //                sumi5[k1] += W[i][k1] * ( 1 / (sigma0[k1] * sigma0[k1]) +  ((resp[i] - mui) / (sigma0[k1] * sigma0[k1] * sigma0[k1])) * (2 * delta[i] - (2 + (resp[i] - mui) / sigma0[k1]) * exp( (resp[i] - mui) / sigma0[k1])  )) ;
    //            }
    //            //sumi5[k1] += sigpennom;
    //            //sumi3[k1] += sigpendenom;
    //            sigma0hat[k1] = sigma0[k1] - pow(0.5, 1) * (1 / sumi5[k1]) * sumi3[k1];
    //        }

    loglike1 = 0.0;
    sat_loglike1 = 0.0;
    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP; k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1] * selection[j][k1];
        mui += alpha0hat[k1];
        deni = pow((1 / sigma0hat[k1]) * exp((resp[i] - mui) / sigma0hat[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0hat[k1])) ;
        //                if (deni < eps2)
        //                    deni = eps2;
        sat_den = pow((1 / sigma0hat[k1]) * exp((resp[i] - 0.0) / sigma0hat[k1]), delta[i]) * exp(- exp((resp[i] - 0.0) / sigma0hat[k1])) ;
        //                if (sat_den < eps2)
        //                    sat_den = eps2;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      loglike1 += log(sumi);
      sat_loglike1 += log(sat_den);
    }

    SUM1 = 0;
    for(k1 = 0; k1 < NCOMP;  k1++)
      for(j = 0; j < NCOV; j++)
        SUM1 += selection[j][k1];

    *loglike = loglike1;

    for (k1 = 0; k1 < NCOMP; k1++)
    {
      for (j = 0;  j < NCOV;  j++)
      {
        mybeta[k1 * NCOV + j] = beta0hat[j][k1];
      }
      myalpha[k1] = alpha0hat[k1];
      mysigma[k1] = sigma0hat[k1];
      mypi[k1] = pi0hat[k1];
    }

    *BIC = loglike1 - 0.5 * SUM1 * log(nsize);
    *EBIC5 = loglike1 - 0.5 * (SUM1) * log(nsize) - 0.5 * (SUM1) * log(NCOV);
    *EBIC1 = loglike1 - 0.5 * (SUM1) * log(nsize) - (SUM1) * log(NCOV);
    *AIC = loglike1 - (SUM1);
    *GCV = (loglike1) / (nsize * pow(1 - SUM1 / nsize, 2));
    *GIC = loglike1 - 0.5 * (SUM1) * log(nsize);
    *MaxEMiter = niter1;
    /*******The E-step of the EM********/

    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP;  k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        deni = pow((1 / sigma0hat[k1])* exp((resp[i] - mui) / sigma0hat[k1]), delta[i]) * exp(-exp((resp[i] - mui) / sigma0hat[k1])) ;
        phi[i][k1] = pi0hat[k1] * deni;
        sumi += phi[i][k1];
      }
      for(k1 = 0; k1 < NCOMP;  k1++){
        W[i][k1] = phi[i][k1] / sumi;
      }
    }

    /**********End of the E-step*******/
    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < nsize; i++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0hat[j][k1];
        mui += alpha0hat[k1];
        predict[k1 * nsize + i] = mui;
        residual[k1 * nsize + i] = resp[i] - mui;
        tau[k1 * nsize + i] = W[i][k1];
      }
    }

  }
#ifdef __cplusplus
}
#endif

/* *************************** Tuning Parameter Weibull ******************* */

#ifdef __cplusplus
extern "C" {
#endif
  void FMR_Weibl_Surv_CwTuneParSel(double *resp,
                                   double *myX,
                                   double *delta,
                                   int *myPenaltyFamily,
                                   double *myridgepen_prop,
                                   int *myNCOMP,
                                   int *myNCOV,
                                   int *mynsize,
                                   double *myinitial_alpha,
                                   double *myinitial_beta,
                                   double *myinitial_sigma,
                                   double *myinitial_pi,
                                   int *myNRmaxiter,
                                   int *myNRportion,
                                   double *myepsconv,
                                   double *gammixportion,
                                   double *optlambda
  )
  {
    int nsize = *mynsize;
    int NCOV = *myNCOV;
    int NCOMP = *myNCOMP;
    int jar = *myPenaltyFamily;

    double ridge1 = *myridgepen_prop;
    double gam = 0.0;
    if(ridge1 != 0.0)
      gam = 1;

    double eps = *myepsconv;
    double GamMP = *gammixportion;

    double NR_maxiter = *myNRmaxiter;
    int NRportion = *myNRportion;
    int alp;


    int i;
    int j;
    int j1;
    int k1;
    int l;
    int check1[1];

    int l1;
    int MaxLim = 80;

    double multX[nsize][NCOV];
    double one_X[nsize][NCOV + 1];

    double optlam[NCOMP];

    double initbeta[NCOV][NCOMP];
    double beta0[NCOV][NCOMP];
    double new_beta0[NCOV][NCOMP];

    double initalpha[NCOMP];
    double alpha0[NCOMP];
    double new_alpha0[NCOMP];

    double initsigma[NCOMP];
    double sigma0[NCOMP];
    double new_sigma0[NCOMP];

    //        double sigpennom = 0.0;
    //        double sigpendenom = 0.0;

    double initpi[NCOMP];
    double pi0[NCOMP];
    double new_pi0[NCOMP];

    double W[nsize][NCOMP];
    double phi[nsize][NCOMP];

    double eps1[NCOMP];

    double vecder[NCOV];
    double vecsigma[NCOV + 1];
    double En[NCOMP];

    double oneXTWY[NCOV + 1];
    double oneXTWX[NCOV + 1][NCOV + 1];
    double selection[NCOV][NCOMP];

    double oneComMat[NCOV + 2][NCOV + 2];
    double onesolution1[NCOV + 2];
    double oneComMatVec[(NCOV + 2) * (NCOV + 2)];

    double sumi;
    double mui;
    double deni;

    double n1[NCOMP];
    double BIC[MaxLim][NCOMP];
    double Max_BIC[NCOMP];
    double lambda1[MaxLim];
    double newloglike1;
    double oldloglike1;
    double jamconvg1;

    double sumi3[NCOMP];
    double sumi5[NCOMP];

    int count1[NCOMP][MaxLim];
    int indx1[NCOMP];

    int niter2;
    char convg2;

    for(l1 = 0; l1 < MaxLim; l1++)
      lambda1[l1] = 0.01 + l1 * 0.01;

    double loglike1;

    double holdveccov[NCOV];
    //    double holdveccom[NCOMP];

    for(k1 = 0; k1 < NCOMP; k1++){
      for(i = 0; i < NCOV; i++){
        new_beta0[i][k1] = beta0[i][k1] = initbeta[i][k1] = myinitial_beta[NCOV * k1 + i];
      }
      new_alpha0[k1] = alpha0[k1] = initalpha[k1] = myinitial_alpha[k1];
      new_sigma0[k1] = sigma0[k1] = initsigma[k1] = myinitial_sigma[k1];
      new_pi0[k1] = pi0[k1] = initpi[k1] = myinitial_pi[k1];
    }

    for (j = 0; j < NCOV; j++){
      for (i = 0; i < nsize; i++){
        multX[i][j] = myX[nsize * j + i];
      }
    }

    for(i = 0; i < nsize; i++){
      one_X[i][0] = 1.0;
      for(j = 1; j < (NCOV + 1); j++){
        one_X[i][j] = multX[i][j - 1];
      }
    }

    for(k1 = 0; k1 < NCOMP; k1++){
      for (i = 0; i < NCOV; i++) {
        holdveccov[i] = initbeta[i][k1];
      }
      if((jar == 1) || (jar == 2))
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (2 * nsize * optlam[k1]);
      else
        eps1[k1] = eps * minimum(holdveccov, NCOV) / (4 * nsize * optlam[k1]);
    }


    for(k1 = 0; k1 < NCOMP; k1++)
      n1[k1] = 0.0;


    for(i = 0; i < nsize; i++){
      sumi = 0.0;
      for(k1 = 0; k1 < NCOMP; k1++){
        mui = 0.0;
        for(j = 0; j < NCOV; j++)
          mui += multX[i][j] * beta0[j][k1];
        mui += alpha0[k1];
        deni = pow((1 / sigma0[k1]) * exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0[k1])) ;
        //                    if (isnan(deni))
        //                        cout << "NAN" << "\t";
        //                    if (deni < eps2)
        //                        deni = eps2;
        phi[i][k1] = pi0[k1] * deni;
        sumi += phi[i][k1];

      }
      for(k1 = 0; k1 < NCOMP; k1++){
        W[i][k1] = phi[i][k1] / sumi;
      }
    }

    //*Start of choosing lambda for each component of the mixture

    for(k1 = 0; k1 < NCOMP; k1++){

      En[k1] = nsize * pow(pi0[k1], GamMP);

      for(l1 = 0; l1 < MaxLim; l1++){

        niter2 = 0;
        convg2 = 'n';
        while((convg2 != 'y') && (niter2 < NR_maxiter)){

          oldloglike1 = 0.0;
          for(i = 0; i < nsize; i++){
            mui = 0.0;
            for(j = 0; j < NCOV; j++)
              mui += multX[i][j] * beta0[j][k1];
            mui += alpha0[k1];
            deni = pow((1 / sigma0[k1]) * exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0[k1])) ;
            //                    if (deni<eps2)
            //                        deni=eps2;
            oldloglike1 += W[i][k1] * log(deni);
          }

          for (j = 0; j < NCOV; j++) {
            holdveccov[j] = beta0[j][k1];
          }

          if(jar == 1)
            for(j = 0; j < NCOV; j++)
              vecder[j] = lambda1[l1];
          else if (jar == 2)
            scadpen(lambda1[l1], vecder, holdveccov, NCOV);
          else if (jar == 3)
            mcppen(lambda1[l1], vecder, holdveccov, NCOV);
          else if (jar == 4)
            sicapen(lambda1[l1], vecder, holdveccov, NCOV);
          else if (jar == 5){
            for(j = 0; j < NCOV; j++)
              vecder[j] = optlam[k1] / (fabs(initbeta[j][k1]) + eps);
          }
          else
            hardpen(optlam[k1], vecder, holdveccov, NCOV);


          for(j = 0; j <(NCOV + 1); j++)
            if(j == 0)
              vecsigma[j] = 0.0;
            else
              vecsigma[j] = vecder[j - 1] / (fabs(beta0[j - 1][k1]) + eps1[k1]);

            /******Constructing the Hessian matrix H ***/
            for(i = 0; i < (NCOV + 1); i++){
              for(j = i; j < (NCOV + 1); j++){
                sumi = 0.0;
                for(l = 0; l < nsize; l++){
                  mui = 0.0;
                  for(j1 = 0; j1 < NCOV; j1++)
                    mui += multX[l][j1] * beta0[j1][k1];
                  mui += alpha0[k1];
                  sumi += - W[l][k1] * one_X[l][i] * one_X[l][j] / (sigma0[k1] * sigma0[k1]) * exp((resp[l] - mui) / sigma0[k1]);
                }
                if(i == j)
                  oneXTWX[i][j] = sumi - En[k1] * vecsigma[j] + ridge1 * log(nsize);
                else
                  oneXTWX[j][i] = oneXTWX[i][j] = sumi;
              }
            }

            oneXTWX[0][0] = oneXTWX[0][0] + (En[k1] * vecsigma[0]) - ridge1 * log(nsize);

            /******Constructing the weigthed vector XTWY***/

            for(j = 0; j < (NCOV + 1); j++){
              sumi = 0.0;
              for(i = 0; i < nsize; i++){
                mui = 0.0;
                for(l = 0; l < NCOV; l++)
                  mui += multX[i][l] * beta0[l][k1];
                mui += alpha0[k1];
                sumi += W[i][k1] * (one_X[i][j] / sigma0[k1]) * (exp((resp[i] - mui) / sigma0[k1]) - delta[i]);
              }
              if (j == 0)
                oneXTWY[j] =  sumi - En[k1] * vecsigma[j];
              else
                oneXTWY[j] =  sumi - En[k1] * vecsigma[j] * beta0[j-1][k1];
            }

            /***In a system Ax=b, adding b to A as its last column**/

            for(i = 0; i < (NCOV + 1); i++)
              for(j = 0; j < (NCOV + 2); j++)
                if(j != (NCOV + 1))
                  oneComMat[i][j] = - oneXTWX[i][j];
                else
                  oneComMat[i][j] = oneXTWY[i];

                for (i = 0; i < (NCOV + 2); i++) {
                  for (j = 0; j < (NCOV + 2); j++) {
                    oneComMatVec[j + i * (NCOV + 2)] = oneComMat[i][j];
                  }
                }

                /**************************************************************/
                /*Solving the system Ax=y to get betahat in the k-th component*/
                /**************************************************************/

                count1[k1][l1] = 0;
                sol( (NCOV + 1) , oneComMatVec, onesolution1, check1);

                alp = 0;
                do {
                  for(j = 0; j < (NCOV + 1); j++)
                    if(j == 0)
                      new_alpha0[k1] = pow(0.5, alp) * onesolution1[j] + alpha0[k1];
                    else
                      new_beta0[j-1][k1] = pow(0.5, alp) * onesolution1[j] + beta0[j-1][k1];

                    sumi5[k1] = 0.0;
                    sumi3[k1] = 0.0;

                    for(i = 0; i < nsize; i++){
                      mui = 0.0;
                      for(l = 0; l < NCOV; l++)
                        mui += multX[i][l] * new_beta0[l][k1];
                      mui += new_alpha0[k1];
                      sumi3[k1] += W[i][k1] * (- delta[i] / sigma0[k1] + ((resp[i] - mui) / (sigma0[k1] * sigma0[k1])) * ( exp( (resp[i] - mui) / sigma0[k1]) - delta[i])   );
                      sumi5[k1] += W[i][k1] * ( delta[i] / (sigma0[k1] * sigma0[k1]) +  ((resp[i] - mui) / (sigma0[k1] * sigma0[k1] * sigma0[k1])) * (2 * delta[i] - (2 + (resp[i] - mui) / sigma0[k1]) * exp( (resp[i] - mui) / sigma0[k1])  )) ;
                    }
                    //sumi5[k1]+=sigpennom;
                    //sumi3[k1]+=sigpendenom;
                    new_sigma0[k1] = sigma0[k1] - pow(0.5, alp) * (1 / sumi5[k1]) * sumi3[k1];
                    //                    new_sigma0[k1] = (new_sigma0[k1] < 0.6)?0.6:new_sigma0[k1];
                    //                    new_sigma0[k1] = (new_sigma0[k1] > 4)?2:new_sigma0[k1];
                    //                    if (k1==0)
                    //                        new_sigma0[k1] =  (new_sigma0[k1]>4)?2:(new_sigma0[k1]);
                    //                    else
                    //                        new_sigma0[k1] =  (new_sigma0[k1]>4)?2:(new_sigma0[k1]);

                    newloglike1 = 0.0;
                    for(i = 0; i < nsize; i++){
                      mui = 0.0;
                      for(j = 0; j < NCOV; j++)
                        mui += multX[i][j] * new_beta0[j][k1];
                      mui += new_alpha0[k1];
                      deni = pow((1 / new_sigma0[k1]) * exp((resp[i] - mui) / new_sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / new_sigma0[k1])) ;
                      //                        if (deni < eps2)
                      //                            deni = eps2;
                      newloglike1 += W[i][k1] * log(deni);
                    }
                    alp++;

                } while ((oldloglike1 > newloglike1) & (alp < NRportion));

                jamconvg1 = 0.0;
                niter2++;

                for(j = 0; j < NCOV; j++)
                  jamconvg1 += pow(new_beta0[j][k1] - beta0[j][k1], 2);
                jamconvg1 += pow(new_pi0[k1] - pi0[k1], 2);
                jamconvg1 += pow(new_alpha0[k1] - alpha0[k1], 2);
                jamconvg1 += pow(new_sigma0[k1] - sigma0[k1], 2);

                convg2 = 'n';
                if (jamconvg1 < eps)
                  convg2 = 'y';

                count1[k1][l1] = 0;
                for(j = 0; j < (NCOV + 1); j++){
                  if(fabs(new_beta0[j - 1][k1]) < 0.2)
                    selection[j - 1][k1] = 0;
                  else
                    selection[j - 1][k1] = 1;
                  count1[k1][l1] += selection[j - 1][k1];
                }

                for(j = 0; j < NCOV; j++)
                  beta0[j][k1] = new_beta0[j][k1];
                alpha0[k1] = new_alpha0[k1];
                sigma0[k1] = new_sigma0[k1];
        }

        loglike1 = 0.0;
        n1[k1] = 0.0;

        for(i = 0; i < nsize; i++){
          mui = 0.0;
          for(j = 0; j < NCOV; j++)
            mui += multX[i][j] * new_beta0[j][k1] * selection[j][k1];
          mui += alpha0[k1];
          deni = pow((1 / sigma0[k1]) * exp((resp[i] - mui) / sigma0[k1]), delta[i]) * exp(- exp((resp[i] - mui) / sigma0[k1])) ;
          loglike1 += W[i][k1] * log(deni);
          n1[k1] += W[i][k1];
        }

        BIC[l1][k1] = loglike1 - 0.5 * (count1[k1][l1]) * log(n1[k1]) - gam * (count1[k1][l1]) * log(NCOV);

        if(l1 == 0){
          Max_BIC[k1] = BIC[l1][k1];
          indx1[k1] = l1;
        }
        else if(BIC[l1][k1] > Max_BIC[k1])
        {
          Max_BIC[k1] = BIC[l1][k1];
          indx1[k1] = l1;
        }
      }//*End of choosing lambda for each component
    }//*End of choosing lambda for both components of the mixture

    for(k1 = 0; k1 < NCOMP; k1++){
      optlambda[k1] = lambda1[indx1[k1]];
    }

  }
#ifdef __cplusplus
}
#endif


