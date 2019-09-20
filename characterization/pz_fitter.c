#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>

#include "MJDSort.h"

/*
 *  fitter routine for pole-zero parameters for HPGe signals
 *    for ORNL MJD analysis codes
 *
 *  returns chisq/DOF
 *  fitted parameters are returned in last 3 parameters of the function call
 */

#define NP 3
#define TRACE_LEN 3000
#define VERB 0

int eval(float *pars, float *sig, int tlo, int thi, float *fit, float deriv[NP][TRACE_LEN]);
int matinv(double *array, int norder, int dim);

/* ======================================================================= */

float pz_fitter(float *sig, int tlo, int thi, int chan, PZinfo *PZI,
                float *lamda1, float *frac2, float *lamda)
{
  double alpha[NP][NP], array[NP][NP], beta[NP];
  float  fit[TRACE_LEN], deriv[NP][TRACE_LEN];
  float  diff, chisq=0.0f, chisq1=0.0f, sigma = 4.0;
  int    j, k, l, n, t;
  int    npars = NP, ndf = thi - tlo - NP;
  float  pars[10] = {0};
  // float errs[NP]={0};


  *lamda1 = *frac2 = *lamda = 0;
  if (thi - tlo < 400) return 0.0f;
  if (thi >= TRACE_LEN) {
    printf("ERROR in pz_fitter: thi = %d is greater than TRACE_LEN = %d\n",
           thi, TRACE_LEN);
    return 0.0f;
  }

  /*
   *  two or three fitted parameters:
   *    pars[0]: signal height
   *    pars[1]: primary decay factor per timestep = timestep/tau
   *    pars[2]: fractional height of overshoot at top of signal
   *    pars[3]: signal height with frac2 fixed at previously determined value
   *    pars[4]: primary decay factor with frac2 fixed at previously determined value
   *  three non-fitted parameters:
   *    pars[5]: resting baseline for this detector, in ADC units
   *    pars[6]: decay of overshoot component = timestep/tau_2
   *    sigma:   rms of signal values to use for chisq
   */

  pars[1] = 0.01 / PZI->tau[chan];    // fractional drop in height over one step
  pars[2] = PZI->frac2[chan];         // fractional height of overshoot at top of signal
  pars[5] = PZI->baseline[chan];      // resting baseline for this detector, in ADC units (fixed)
  pars[6] = 0.01 / PZI->tau2[chan];   // decay per timestep of overshoot component (fixed)
  sigma   = PZI->bl_rms[chan];        // mean RMS of signal baseline, used to determine chisq
  if (sigma < 0.1) sigma = 4.0;

  /* evaluate fit, alpha & beta matrices, & chisq */
  eval(pars, sig, tlo, thi, fit, deriv);
  pars[3] = pars[0];  // save for second fit
  pars[4] = pars[1];  // save for second fit
  ndf = thi - tlo - NP;
  for (j = 0; j < npars; ++j) {
    beta[j] = 0.0f;
    for (k = 0; k <= j; ++k) {
      alpha[k][j] = 0.0f;
    }
  }

  for (t = tlo; t < thi; t++) {
    diff = sig[t] - fit[t];
    chisq1 += diff * diff;
    for (l=0; l<npars; l++) {
      beta[l] += diff * deriv[l][t];
      for (n=0; n<=l; n++) {
        alpha[n][l] += (double)deriv[l][t] * (double)deriv[n][t];
      }
    }
  }

  chisq1 /= (float) ndf * sigma*sigma;
  /* invert modified curvature matrix to find new parameters */
  for (j = 0; j < npars; ++j) {
    if (alpha[j][j] * alpha[j][j] == 0.) {
      printf("Cannot fit PZ pars - diagonal element no. %d equal to zero.\n", j);
      return -1.0f;
    }
  }
  array[0][0] = 1.0f;
  for (j = 1; j < npars; ++j) {
    for (k = 0; k < j; ++k) {
      array[k][j] = alpha[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
      array[j][k] = array[k][j];
    }
    array[j][j] = 1.0f;
  }
  matinv(array[0], npars, NP);
  /* calculate new par. values and chisq */
  for (j = 0; j < npars; ++j) {
    for (k = 0; k < npars; ++k) {
      pars[j] += beta[k] * array[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
    }
  }
  /*  calculate errors??
    for (j = 0; j < npars; ++j) {
    if (array[j][j] < 0.) array[j][j] = 0.;
    errs[j] = sqrt(chisq * array[j][j] / alpha[j][j]);
    if (VERB) printf ("%d %f %f\n", j, pars[j], errs[j]); fflush(stdout);
    }
  */
  if (pars[1] < 0 || pars[2] > 0.1) {
    return -1.0f;
  }
  *lamda1 = pars[1];
  *frac2  = pars[2];

  /* calculate new chisq */
  eval(pars, sig, tlo, thi, fit, deriv);
  for (t = tlo; t < thi; t++) {
    diff = sig[t] - fit[t];
    chisq += diff * diff;
  }
  chisq /= (float) ndf * sigma*sigma;

  if (VERB) {
    printf("*** chisq, chisq1 = %.4f %.4f  ***  pars = %f %f %f\n",
           chisq, chisq1, pars[0], pars[1]*1000.0, pars[2]);
    fflush(stdout);
  }

  /* now redo fit with only two parameters (i.e. fix frac2 at original value) */
  array[0][0] = array[1][1] = 1.0f;
  array[0][1] = array[1][0] = alpha[0][1] / sqrt(alpha[0][0] * alpha[1][1]);
  matinv(array[0], 2, NP);
  /* calculate new par. values and chisq */
  for (j = 0; j < 2; ++j) {
    for (k = 0; k < 2; ++k) {
      pars[j+NP] += beta[k] * array[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
    }
  }
  /*  calculate errors??
    for (j = 0; j < npars; ++j) {
    if (array[j][j] < 0.) array[j][j] = 0.;
    errs[j] = sqrt(chisq * array[j][j] / alpha[j][j]);
    if (VERB) printf ("%d %f %f\n", j, pars[j], errs[j]); fflush(stdout);
    }
  */
  if (0) {  // chisq calculation with fixed frac2; not needed?
    chisq = 0;
    pars[0] = pars[3];   // newly fitted values
    pars[1] = pars[4];   // newly fitted values
    eval(pars, sig, tlo, thi, fit, deriv);
    for (t = tlo; t < thi; t++) {
      diff = sig[t] - fit[t];
      chisq += diff * diff;
    }
    chisq /= (float) (ndf+1) * sigma*sigma;  // +1 since one fixed parameter
  }

  *lamda = pars[4];
  return chisq;
} /* fitter */

/* ======================================================================= */

/*
 *  two or three fitted parameters:
 *    pars[0]: signal height
 *    pars[1]: primary decay factor per time step = timestep/tau
 *    pars[2]: fractional height of overshoot at top of signal
 *    pars[3]: signal height with frac2 fixed at previouslt determined value
 *    pars[4]: primary decay factor with frac2 fixed at previouslt determined value
 *  three non-fitted parameters:
 *    pars[5]: resting baseline for this detector, in ADC uints
 *    pars[6]: decay of overshoot component = timestep/tau_2
 *    pars[7]: final chisq/DOF
 *   [pars[8]: final chisq/DOF with frac2 fixed]
 *
 */

int eval(float *pars, float *sig, int tlo, int thi, float *fit, float deriv[NP][TRACE_LEN])
{

  double decay = pars[1], decay2 = pars[6];   // fractional drop in height over one time step
  float  baseline = pars[5], frac2 = pars[2];
  float  base=0, e2, e3=0, e4;
  int    i;
  // float  fsig[TRACE_LEN];

  /* find starting baseline */
  for (i=20; i<40; i++) base += sig[i];
  base /= 20.0;
  if (base > baseline) base = baseline;

  deriv[1][tlo-1] = 0;

  /* start at beginning to do PZ correction, going as far as tlo */
  e2 = sig[30];
  for (i = 31; i < tlo; i++) {
    e3 += sig[i] - e2 - e3*decay2;
    e2  = sig[i];
  }

  float s = (sig[tlo-1] + sig[tlo] + sig[tlo+1]) / 3.0;
  pars[0] = s - base - e3 * frac2; // initial esitimate of signal height
  if (0) printf("\n Fit from %d to %d, pars %f %f %f %f %f  base %f\n\n",
                tlo, thi, pars[0], pars[1], pars[2], pars[5], pars[6], base);
  e4 = pars[0];

  /* continue to do PZ correction, now calculating derivatives */
  for (; i <= thi; i++) {
    e3 += sig[i] - e2 - e3*decay2;
    e2  = sig[i];
    deriv[0][i] = e4/pars[0];
    deriv[1][i] = -e4*(double)(i-tlo+1);
    deriv[2][i] = e3;
    fit[i] = base + e4 + frac2*e3;
    e4 *= (1.0 - decay);
  }

  return 0;
} /* eval */


/* ======================================================================= */
int matinv(double *array, int norder, int dim)
{
  double amax, save;
  int i, j, k, ik[4000], jk[4000];

  for (k = 0; k < norder; ++k) {
    /* find largest element array(i,j) in rest of matrix */
    amax = 0.f;
    while (1) {
      for (i = k; i < norder; ++i) {
	for (j = k; j < norder; ++j) {
	  if (fabs(amax) < fabs(array[i + j*dim])) {
	    amax = array[i + j*dim];
	    ik[k] = i;
	    jk[k] = j;
	  }
	}
      }
      if (amax == 0.f) {
	printf("\nmatinv: k = %d, npars = %d!\n", k, norder);
	return 0;
      }
      /* interchange rows and columns to put amax in array(k,k) */
      i = ik[k];
      if (i < k) continue;
      if (i > k) {
	for (j = 0; j < norder; ++j) {
	  save = array[k + j*dim];
	  array[k + j*dim] = array[i + j*dim];
	  array[i + j*dim] = -save;
	}
      }
      j = jk[k];
      if (j >= k) break;
    }
    if (j > k) {
      for (i = 0; i < norder; ++i) {
	save = array[i + k*dim];
	array[i + k*dim] = array[i + j*dim];
	array[i + j*dim] = -save;
      }
    }
    /* accumulate elements of inverse matrix */
    for (i = 0; i < norder; ++i) {
      if (i != k) array[i + k*dim] = -array[i + k*dim] / amax;
    }
    for (i = 0; i < norder; ++i) {
      for (j = 0; j < norder; ++j) {
	if (i != k && j != k)
	  array[i + j*dim] += array[i + k*dim] * array[k + j*dim];
      }
    }
    for (j = 0; j < norder; ++j) {
      if (j != k) array[k + j*dim] /= amax;
    }
    array[k + k*dim] = 1.f / amax;
  }
  /* restore ordering of matrix */
  for (k = norder-1; k >= 0; --k) {
    j = ik[k];
    if (j > k) {
      for (i = 0; i < norder; ++i) {
	save = array[i + k*dim];
	array[i + k*dim] = -array[i + j*dim];
	array[i + j*dim] = save;
      }
    }
    i = jk[k];
    if (i > k) {
      for (j = 0; j < norder; ++j) {
	save = array[k + j*dim];
	array[k + j*dim] = -array[i + j*dim];
	array[i + j*dim] = save;
      }
    }
  }
  return 0;
} /* matinv */
