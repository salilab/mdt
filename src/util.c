/** \file util.c           MDT utility functions
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <math.h>
#include "modeller.h"
#include "mod_dynmem.h"
#include "util.h"

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
int *mdt_start_indices(const struct mdt_type *mdt)
{
  int i, *indf;

  indf = dmalloc(sizeof(int) * mdt->nfeat);
  for (i = 0; i < mdt->nfeat; i++) {
    indf[i] = f_int1_get(&mdt->istart, i);
  }
  return indf;
}

/** Calculate the weights in the smoothing procedure for combining the
    a priori pdf with the experimental pdf. */
void weights(float weight, int nbins, float norm, float *w1, float *w2)
{
  if (weight > 1e-10) {
    *w1 = 1. / (1. + norm / (weight * (float)nbins));
  } else {
    *w1 = 0.;
  }
  *w2 = 1. - *w1;
}


/** Return the index in the MDT pdf of the point whose feature bin
    indices are indf. */
int indmdt(const int *indf, const struct mdt_type *mdt)
{
  int i, ind;
  ind = indf[mdt->nfeat - 1] - f_int1_get(&mdt->istart, mdt->nfeat - 1);

  for (i = mdt->nfeat - 2; i >= 0; i--) {
    int indval = indf[i] - f_int1_get(&mdt->istart, i);
    ind += f_int1_get(&mdt->stride, i) * indval;
  }
  return ind;
}


/** Update the indices for the next point in the MDT. Return false if no
    more points are available. */
int roll_ind(int indf[], const int istart[], const int iend[], int nfeat)
{
  int i = nfeat - 1;
  while (i >= 0) {
    if (indf[i] + 1 <= iend[i]) {
      indf[i]++;
      return 1;
    } else if (i == 0) {
      return 0;
    } else {
      indf[i] = istart[i];
      i--;
    }
  }
  return 0;
}

/** Get the number of bins in the 1 or 2 dependent features */
void get_binx_biny(int dimensions, const struct mdt_type *mdt,
                   const char *routine, int *nbinx, int *nbiny, int *ierr)
{
  if (dimensions < 1 || dimensions > 2 || dimensions > mdt->nfeat) {
    modlogerror(routine, ME_VALUE,
                "'dimensions' is %d; it must be either 1 or 2, and not more "
                "than the dimensionality of this MDT (%d)", dimensions,
                mdt->nfeat);
    *ierr = 1;
  } else {
    *ierr = 0;
    if (dimensions == 1) {
      *nbinx = f_int1_get(&mdt->nbins, mdt->nfeat - 1);
      *nbiny = 1;
    } else {
      *nbinx = f_int1_get(&mdt->nbins, mdt->nfeat - 1);
      *nbiny = f_int1_get(&mdt->nbins, mdt->nfeat - 2);
    }
  }
}


/** Return the sum of a set. */
double get_sum(const double bin[], int nbins)
{
  int i;
  double sum = 0.;
  for (i = 0; i < nbins; i++) {
    sum += bin[i];
  }
  return sum;
}


/** Return the entropy of a set. */
double entrp1(const double frq[], int nbinx)
{
  int i;
  double ent;

  ent = 0.;
  for (i = 0; i < nbinx; i++) {
    if (frq[i] > 1.0e-37) {
      ent -= frq[i] * log(frq[i]);
    }
  }
  return ent;
}
