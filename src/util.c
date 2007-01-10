/** \file util.c           MDT utility functions
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <math.h>
#include <assert.h>
#include "modeller.h"
#include "mod_dynmem.h"
#include "util.h"
#include "num_recipes.h"

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


/** Roll n indices in ind so that all combinations of n different indices
    are generated, where the possible values for each index are 1+x to nmax-y.

    For example, for n=3:

    for (i = 0, i < nmax-2; i++) {
      for (j = i+1, j < nmax-1; j++) {
        for (k = j+1, k < nmax; k++) {

    *ind should be NULL on the first call to this routine, and it will be
    initialized (the user should free it when finished).
    Returns false if no more indices are available. */
int roll_ind_comb(int **ind, int n, int nmax)
{
  int i;
  int *indr;

  if (n == 0) {
    return 0;
  } else if (*ind == NULL) {
    *ind = indr = dmalloc(sizeof(int) * n);
    for (i = 0; i < n; i++) {
      indr[i] = i;
    }
    return 1;
  } else {
    indr = *ind;
    for (i = n; i > 0; i--) {
      if (indr[i] < nmax - (n-i)) {
        int k;
        indr[i]++;
        for (k = i + 1; k < n; k++) {
          indr[k] = indr[k-1] + 1;
        }
        return 1;
      }
    }
    return 0;
  }
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

/** Allocate and set control arrays for a subset of the MDT features */
static void setup_mdt_feature_arrays(const struct mdt_type *mdt,
                                     const int ifeat[], int nfeat, int **istart,
                                     int **ival, int **nbins)
{
  int i;
  *istart = dmalloc(sizeof(int) * nfeat);
  *ival = dmalloc(sizeof(int) * nfeat);
  *nbins = dmalloc(sizeof(int) * nfeat);
  for (i = 0; i < nfeat; i++) {
    (*nbins)[i] = f_int1_get(&mdt->nbins, ifeat[i]);
    (*ival)[i] = (*istart)[i] = f_int1_get(&mdt->istart, ifeat[i]);
  }
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
  setup_mdt_feature_arrays(mdt, i_feat_var, n_feat_var, &istart, &i_val_var,
                           &var_feature_bins);

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

/** Return entropy of p(x/y,z,...) where y,z are the independent features.
    See pages 480-483 in Numerical Recipes for equations. */
double entrp2(double summdt, const int i_feat_fix[], const struct mdt_type *mdt,
              int n_feat_fix, int nbinx, float sumi[])
{
  static const char *routine = "entrp2";
  static const double small = 1.0e-3, tiny = 1.0e-6;
  int i, *fix_feature_bins, *i_val_fix, *istart;
  double spjk, entrp, sumfrq, *fijk;

  /* initialize the column and row totals; only for the chi^2 calc later on: */
  for (i = 0; i < nbinx; i++) {
    sumi[i] = 0.;
  }

  fijk = dmalloc(sizeof(double) * nbinx);
  /* set the control arrays for fixed (and non-independent) feature values: */
  setup_mdt_feature_arrays(mdt, i_feat_fix, n_feat_fix, &istart, &i_val_fix,
                           &fix_feature_bins);

  /* for each possible combination of values of the current combination
     of independent features (all j,k,...) */
  spjk = entrp = 0.;
  do {
    /* get pdf p(x/y,z,...) for the current combination of values of the
       current combination of independent features */
    getfrq(mdt, i_feat_fix, n_feat_fix, i_val_fix, nbinx, fijk);

    /* over all values of the dependent variable */
    sumfrq = get_sum(fijk, nbinx);

    /* only for the chi^2 calculation later on: */
    for (i = 0; i < nbinx; i++) {
      sumi[i] += fijk[i];
    }

    if (sumfrq > tiny) {
      double e2, pjk;
      int ix;
      /* for all values of the dependent variable */
      e2 = 0.;
      for (ix = 0; ix < nbinx; ix++) {
        double p = fijk[ix] / sumfrq;
        /* note that: lim_{p->0} p*ln(p) = 0 */
        if (p > tiny) {
          e2 += p * log(p);
        }
      }

      /* sumfrq/summdt is a weight for this contribution to total entropy:
         it equals to the probability of occurence of this particular
         combination of independent feature values in the sample: */
      pjk = sumfrq / summdt;
      entrp -= pjk * e2;
      spjk += pjk;
    }

  } while (roll_ind(i_val_fix, istart, fix_feature_bins, n_feat_fix));

  if (fabs(spjk - 1.0) > small) {
    modlogwarning(routine, "Sum of p(ijk) (%.8f) is not equal to 1.0", spjk);
  }

  /* clean up */
  free(fijk);
  free(fix_feature_bins);
  free(i_val_fix);
  free(istart);

  return entrp;
}

/** Get the chi^2, etc for pdf p(x/y,z,...) */
double chisqr(double summdt, const int i_feat_fix[], const struct mdt_type *mdt,
              int n_feat_fix, int nbinx, float sumi[], double *df, double *prob,
              double *ccc, double *cramrv, int *ierr)
{
  static const char *routine = "chiqrt";
  static const double small = 1.0e-6, tiny = 1.0e-30, fract = 0.001;
  int i, *fix_feature_bins, *i_val_fix, *istart;
  float s, sumjmdt;
  double *fijk, chi2;
  int nni, nnj, minnn;

  *ierr = 0;

  /* a self-consistency test and getting nni: */
  nni = 0;
  s = 0.;
  for (i = 0; i < nbinx; i++) {
    s += sumi[i];
    if (sumi[i] > small) {
      nni++;
    }
  }
  if (fabs((s - summdt) / (s + tiny)) > fract) {
    modlogerror(routine, ME_GENERIC,
                "Inconsistency in chi-square calculation: %f %f", s, summdt);
    *ierr = 1;
    return 0.;
  }

  fijk = dmalloc(sizeof(double) * nbinx);
  /* set the control arrays for non-independent feature values: */
  setup_mdt_feature_arrays(mdt, i_feat_fix, n_feat_fix, &istart, &i_val_fix,
                           &fix_feature_bins);

  /* do the chi^2 sum: */
  nnj = 0;
  sumjmdt = 0.;
  chi2 = 0.;
  do {
    double sumj;

    /* get pdf p(x/y,z,...) for the current combination of values of the
       current combination of independent features */
    getfrq(mdt, i_feat_fix, n_feat_fix, i_val_fix, nbinx, fijk);

    sumj = get_sum(fijk, nbinx);
    sumjmdt += sumj;
    if (sumj > small) {
      nnj++;
    }

    for (i = 0; i < nbinx; i++) {
      double expctd, actual;
      expctd = sumj * sumi[i] / summdt;
      actual = fijk[i];
      chi2 += (actual - expctd) * (actual - expctd) / (expctd + tiny);
    }
  } while (roll_ind(i_val_fix, istart, fix_feature_bins, n_feat_fix));

  /* clean up */
  free(fijk);
  free(fix_feature_bins);
  free(i_val_fix);
  free(istart);

  /* test for the sum over j: */
  if (fabs((sumjmdt - summdt) / (summdt + tiny)) > fract) {
    modlogerror(routine, ME_GENERIC,
                "Inconsistency in chi-square calculation: %f %f", sumjmdt,
                summdt);
    *ierr = 1;
    return 0.;
  }

  /* degrees of freedom for chi^2: */
  *df = nni * nnj - nni - nnj + 1;

  /* significance level of the chi^2 (small values indicating significant
     association): */
  if (*df > small) {
    *prob = gammq(0.5 * (*df), 0.5 * chi2);
  } else {
    modlogwarning(routine, "There are zero degrees of freedom; "
                  "chi-square calculation is meaningless.");
    *prob = -1.0;
  }

  /* Cramer V (or phi statistic) for the strength of the association: */
  minnn = (nni < nnj ? nni : nnj);
  *cramrv = sqrt(chi2 / (summdt * (minnn - 1)));

  /* contingency coefficient C for the strength of the association: */
  *ccc = sqrt(chi2 / (chi2 + summdt));

  return chi2;
}
