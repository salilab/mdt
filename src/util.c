/** \file util.c           MDT utility functions
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <math.h>
#include <assert.h>
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

/** Given a set of independent (fixed) variables, find the complementary set
    of dependent features, excluding the last. Return the indices and the number
    of the dependent (variable) features.
   */
static int *complem(const int i_feat_fix[], int n_feat_fix, int numb_features,
                    int *n_feat_var)
{
  int i, *i_feat_var;

  i_feat_var = dmalloc(sizeof(int) * numb_features);
  *n_feat_var = 0;
  for (i = 0; i < numb_features - 1; i++) {
    int j, match = 0;
    for (j = 0; j < n_feat_fix; j++) {
      if (i_feat_fix[j] == i) {
        match = 1;
        break;
      }
    }
    if (!match) {
      i_feat_var[*n_feat_var] = i;
      (*n_feat_var)++;
    }
  }
  assert(*n_feat_var + n_feat_fix + 1 == numb_features);
  return i_feat_var;
}


/** Calculate the pdf p(x/independent features) summed over all
    independent features and their values except for the n_feat_fix fixed
    independent features. */
void getfrq(const struct mdt_type *mdt, const int i_feat_fix[], int n_feat_fix,
            const int i_val_fix[], int nbinx, double frq[])
{
  int *indf, *i_val_var, *i_feat_var, *var_feature_bins, *istart, n_feat_var,
      i, i1;
  double *bin;

  /* Find the dependent features */
  i_feat_var = complem(i_feat_fix, n_feat_fix, mdt->nfeat, &n_feat_var);

  indf = mdt_start_indices(mdt);

  /* set the fixed feature values ordered for indmdt(): */
  for (i = 0; i < n_feat_fix; i++) {
    indf[i_feat_fix[i]] = i_val_fix[i];
  }

  /* set the control arrays for non-independent feature values: */
  istart = dmalloc(sizeof(int) * n_feat_var);
  i_val_var = dmalloc(sizeof(int) * n_feat_var);
  var_feature_bins = dmalloc(sizeof(int) * n_feat_var);
  for (i = 0; i < n_feat_var; i++) {
    var_feature_bins[i] = f_int1_get(&mdt->nbins, i_feat_var[i]);
    i_val_var[i] = istart[i] = f_int1_get(&mdt->istart, i);
  }

  /* clear the output frq array: */
  for (i = 0; i < nbinx; i++) {
    frq[i] = 0.;
  }

  bin = f_double1_pt(&mdt->bin);
  /* for all combinations of the values of the variable features: */
  do {
    /* copy new values of variable features into the indf array
       (the fixed values are set above) */
    for (i = 0; i < n_feat_var; i++) {
      indf[i_feat_var[i]] = i_val_var[i];
    }

    /* for each x, sum mdt over all possible values for n_feat_var features */
    i1 = indmdt(indf, mdt);
    for (i = 0; i < nbinx; i++) {
      frq[i] += bin[i1 + i];
    }
  } while (roll_ind(i_val_var, istart, var_feature_bins, n_feat_var));

  /* Clean up */
  free(i_feat_var);
  free(istart);
  free(i_val_var);
  free(indf);
  free(var_feature_bins);
}
