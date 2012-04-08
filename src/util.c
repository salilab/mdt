/** \file util.c           MDT utility functions
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#include <math.h>
#include <assert.h>
#include <glib.h>
#include "modeller.h"
#include "util.h"
#include "num_recipes.h"

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
int *mdt_start_indices(const struct mod_mdt *mdt)
{
  int i, *indf;

  indf = g_malloc(sizeof(int) * mdt->nfeat);
  for (i = 0; i < mdt->nfeat; i++) {
    indf[i] = mdt->features[i].istart;
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
int indmdt(const int *indf, const struct mod_mdt *mdt)
{
  int i, ind, nfeat = mdt->nfeat;
  const struct mod_mdt_feature *feat = mdt->features;
  ind = indf[nfeat - 1] - feat[nfeat - 1].istart;

  for (i = nfeat - 2; i >= 0; i--) {
    int indval = indf[i] - feat[i].istart;
    ind += feat[i].stride * indval;
  }
  return ind;
}


/** Return the index in the MDT pdf of the point whose feature bin
    indices are indf, using the stride and istart arrays. */
int indmdt_full(const int *indf, const int stride[], int nfeat,
                const int istart[])
{
  int i, ind;
  ind = indf[nfeat - 1] - istart[nfeat - 1];

  for (i = nfeat - 2; i >= 0; i--) {
    int indval = indf[i] - istart[i];
    ind += stride[i] * indval;
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
    } else {
      indf[i] = istart[i];
      i--;
    }
  }
  return 0;
}

/** Update the indices for the next point in the MDT. Return false if no
    more points are available. */
int roll_ind_mdt(int indf[], const struct mod_mdt *mdt, int nfeat)
{
  int i = nfeat - 1;
  while (i >= 0) {
    if (indf[i] + 1 <= mdt->features[i].iend) {
      indf[i]++;
      return 1;
    } else if (i == 0) {
      return 0;
    } else {
      indf[i] = mdt->features[i].istart;
      i--;
    }
  }
  return 0;
}

/** Like roll_ind(), but only for the selected inds[n_inds] features */
int roll_inds(int indf[], const struct mod_mdt *mdt, const int inds[],
              int n_inds)
{
  int iind, i;
  assert(n_inds > 0 && n_inds <= mdt->nfeat);
  iind = n_inds - 1;
  while (iind >= 0) {
    i = inds[iind];
    assert(i >= 0 && i < mdt->nfeat);
    if (indf[i] + 1 <= mdt->features[i].iend) {
      indf[i]++;
      return 1;
    } else if (i == 0) {
      return 0;
    } else {
      indf[i] = mdt->features[i].istart;
      iind--;
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

  if (*ind == NULL) {
    /* In the case where n == 0, allocate a dummy size 1 array so that the
       pointer is not NULL */
    *ind = indr = g_malloc(sizeof(int) * (n == 0 ? 1 : n));
    for (i = 0; i < n; i++) {
      indr[i] = i;
    }
    return 1;
  } else if (n == 0) {
    return 0;
  } else {
    indr = *ind;
    for (i = n - 1; i >= 0; i--) {
      if (indr[i] < nmax - (n - i)) {
        int k;
        indr[i]++;
        for (k = i + 1; k < n; k++) {
          indr[k] = indr[k - 1] + 1;
        }
        return 1;
      }
    }
    return 0;
  }
}


/** Get the number of bins in the 1 or 2 dependent features. Return TRUE on
    success. */
gboolean get_binx_biny(int dimensions, const struct mod_mdt *mdt,
                       const char *routine, int *nbinx, int *nbiny,
                       GError **err)
{
  if (dimensions < 1 || dimensions > 2 || dimensions > mdt->nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: 'dimensions' is %d; it must be either 1 or 2, and not "
                "more than the dimensionality of this MDT (%d)", routine,
                dimensions, mdt->nfeat);
    return FALSE;
  } else {
    if (dimensions == 1) {
      *nbinx = mdt->features[mdt->nfeat - 1].nbins;
      *nbiny = 1;
    } else {
      *nbinx = mdt->features[mdt->nfeat - 1].nbins;
      *nbiny = mdt->features[mdt->nfeat - 2].nbins;
    }
    return TRUE;
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

/** Return the sum of an MDT. */
double get_mdt_sum(const struct mod_mdt *mdt)
{
  int i;
  double sum = 0.;
  for (i = 0; i < mdt->nelems; ++i) {
    sum += mod_mdt_bin_get(mdt, i);
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
    of dependent features, excluding the last dimensions features. Return the
    indices and the number of the dependent (variable) features.
   */
static int *complem(const int i_feat_fix[], int n_feat_fix, int numb_features,
                    int dimensions, int *n_feat_var)
{
  int i, *i_feat_var;

  i_feat_var = g_malloc(sizeof(int) * numb_features);
  *n_feat_var = 0;
  for (i = 0; i < numb_features - dimensions; i++) {
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
  assert(*n_feat_var + n_feat_fix + dimensions == numb_features);
  return i_feat_var;
}

/** Allocate and set control arrays for a subset of the MDT features */
static void setup_mdt_feature_arrays(const struct mod_mdt *mdt,
                                     const int ifeat[], int nfeat,
                                     int **istart, int **ival, int **nbins)
{
  int i;
  *istart = g_malloc(sizeof(int) * nfeat);
  *ival = g_malloc(sizeof(int) * nfeat);
  *nbins = g_malloc(sizeof(int) * nfeat);
  for (i = 0; i < nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->features[ifeat[i]];
    (*nbins)[i] = feat->nbins;
    (*ival)[i] = (*istart)[i] = feat->istart;
  }
}

/** Calculate the pdf p(x/independent features) summed over all
    independent features and their values except for the n_feat_fix fixed
    independent features. */
void getfrq(const struct mod_mdt *mdt, const int i_feat_fix[], int n_feat_fix,
            const int i_val_fix[], int dimensions, int nbinx, double frq[])
{
  int *indf, *i_val_var, *i_feat_var, *var_feature_bins, *istart, n_feat_var,
      i, i1;

  /* Find the dependent features */
  i_feat_var = complem(i_feat_fix, n_feat_fix, mdt->nfeat, dimensions,
                       &n_feat_var);

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

  /* for all combinations of the values of the variable features: */
  do {
    /* copy new values of variable features into the indf array
       (the fixed values are set above) */
    for (i = 0; i < n_feat_var; i++) {
      indf[i_feat_var[i]] = i_val_var[i];
    }

    /* for each x, sum mdt over all possible values for n_feat_var features */
    i1 = indmdt(indf, mdt);
    assert(i1 + nbinx - 1 < mdt->nelems);
    for (i = 0; i < nbinx; i++) {
      frq[i] += mod_mdt_bin_get(mdt, i1 + i);
    }
  } while (roll_ind(i_val_var, istart, var_feature_bins, n_feat_var));

  /* Clean up */
  g_free(i_feat_var);
  g_free(istart);
  g_free(i_val_var);
  g_free(indf);
  g_free(var_feature_bins);
}

/** Return entropy of p(x/y,z,...) where y,z are the independent features.
    See pages 480-483 in Numerical Recipes for equations. */
double entrp2(double summdt, const int i_feat_fix[],
              const struct mod_mdt *mdt, int n_feat_fix, int nbinx,
              float sumi[])
{
  static const char *routine = "entrp2";
  static const double small = 1.0e-3, tiny = 1.0e-6;
  int i, *fix_feature_bins, *i_val_fix, *istart;
  double spjk, entrp, sumfrq, *fijk;

  /* initialize the column and row totals; only for the chi^2 calc later on: */
  for (i = 0; i < nbinx; i++) {
    sumi[i] = 0.;
  }

  fijk = g_malloc(sizeof(double) * nbinx);
  /* set the control arrays for fixed (and non-independent) feature values: */
  setup_mdt_feature_arrays(mdt, i_feat_fix, n_feat_fix, &istart, &i_val_fix,
                           &fix_feature_bins);

  /* for each possible combination of values of the current combination
     of independent features (all j,k,...) */
  spjk = entrp = 0.;
  do {
    /* get pdf p(x/y,z,...) for the current combination of values of the
       current combination of independent features */
    getfrq(mdt, i_feat_fix, n_feat_fix, i_val_fix, 1, nbinx, fijk);

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
    mod_logwarning(routine, "Sum of p(ijk) (%.8f) is not equal to 1.0", spjk);
  }

  /* clean up */
  g_free(fijk);
  g_free(fix_feature_bins);
  g_free(i_val_fix);
  g_free(istart);

  return entrp;
}

/** Get the chi^2, etc for pdf p(x/y,z,...) */
double chisqr(double summdt, const int i_feat_fix[],
              const struct mod_mdt *mdt, int n_feat_fix, int nbinx,
              float sumi[], double *df, double *prob, double *ccc,
              double *cramrv, GError **err)
{
  static const char *routine = "chiqrt";
  static const double small = 1.0e-6, tiny = 1.0e-30, fract = 0.001;
  int i, *fix_feature_bins, *i_val_fix, *istart;
  float s, sumjmdt;
  double *fijk, chi2;
  int nni, nnj, minnn;

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
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "%s: Inconsistency in chi-square calculation: %f %f", routine,
                s, summdt);
    return 0.;
  }

  fijk = g_malloc(sizeof(double) * nbinx);
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
    getfrq(mdt, i_feat_fix, n_feat_fix, i_val_fix, 1, nbinx, fijk);

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
  g_free(fijk);
  g_free(fix_feature_bins);
  g_free(i_val_fix);
  g_free(istart);

  /* test for the sum over j: */
  if (fabs((sumjmdt - summdt) / (summdt + tiny)) > fract) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "%s: Inconsistency in chi-square calculation: %f %f", routine,
                sumjmdt, summdt);
    return 0.;
  }

  /* degrees of freedom for chi^2: */
  *df = nni * nnj - nni - nnj + 1;

  /* significance level of the chi^2 (small values indicating significant
     association): */
  if (*df > small) {
    *prob = gammq(0.5 * (*df), 0.5 * chi2);
  } else {
    mod_logwarning(routine, "There are zero degrees of freedom; "
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

/** Make the stride array for faster indmdt lookup */
void make_mdt_stride(struct mod_mdt *mdt)
{
  int i, nelems;

  mdt->features[mdt->nfeat - 1].stride = 1;
  for (i = mdt->nfeat - 2; i >= 0; i--) {
    mdt->features[i].stride = mdt->features[i + 1].stride
        * mdt->features[i + 1].nbins;
  }

  /* number of elements in the full MDT array: */
  nelems = mdt->features[0].stride * mdt->features[0].nbins;
  assert(mdt->nelems == nelems);
}

/** Make the stride array from the nbins array, and return the size of
    the MDT. */
int make_mdt_stride_full(const int nbins[], int nfeat, int stride[])
{
  int i;

  stride[nfeat - 1] = 1;
  for (i = nfeat - 2; i >= 0; i--) {
    stride[i] = stride[i + 1] * nbins[i + 1];
  }
  return stride[0] * nbins[0];
}

/** Open a file, and uncompress it if necessary. */
struct mod_file *mdt_open_file(const char *path, const char *mode, GError **err)
{
  struct mod_file *fh = mod_file_open(path, mode);
  if (!fh) {
    GError *moderr = mod_error_get();
    if (moderr) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_IO, "%s", moderr->message);
      g_error_free(moderr);
    }
  }
  return fh;
}

/** Close an open file, and do any other necessary tidy-up if it was
    compressed. The initial value of err is used (if an error was already
    set, it is not modified, but emergency cleanup is done here). */
gboolean mdt_close_file(struct mod_file *fh, GError **err)
{
  int ierr;
  if (err && *err) {
    ierr = 1;
  } else {
    ierr = 0;
  }
  mod_error_clear();
  mod_file_close(fh, &ierr);
  if (ierr) {
    GError *moderr = mod_error_get();
    if (moderr) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_IO, "%s", moderr->message);
      g_error_free(moderr);
    }
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Convert a Modeller error into a GError */
void handle_modeller_error(GError **err)
{
  GError *moderr = mod_error_get();
  g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "%s", moderr->message);
  g_error_free(moderr);
}

/** Check that a given feature type is within range. Return FALSE and set the
    error indicator if it is not. */
gboolean check_feature_type(int ifeat, const struct mdt_library *mlib,
                            GError **err)
{
  if (ifeat < 1 || ifeat > mlib->base.nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                "Feature type %d out of range 1-%d", ifeat, mlib->base.nfeat);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Get the position in the MDT bin array of the given set of indices. Return
    TRUE on success. */
gboolean get_bin_index(const struct mod_mdt *mdt, const int indices[],
                       int n_indices, int *bin_index, GError **err)
{
  int i, *indf, indx;
  if (n_indices != mdt->nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "number of indices (%d) must match dimension of MDT (%d)",
                n_indices, mdt->nfeat);
    return FALSE;
  }
  indf = g_malloc(sizeof(int) * n_indices);
  for (i = 0; i < n_indices; i++) {
    /* count negative indices from the end of the feature, as in Python */
    if (indices[i] < 0) {
      indf[i] = mdt->features[i].iend + 1 + indices[i];
    } else {
      indf[i] = indices[i] + 1;
    }
    if (indf[i] < mdt->features[i].istart || indf[i] > mdt->features[i].iend) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                  "index (%d) out of range (%d<=index<=%d) in dimension %d",
                  indf[i] - 1, mdt->features[i].istart - 1,
                  mdt->features[i].iend - 1, i);
      g_free(indf);
      return FALSE;
    }
  }
  indx = indmdt(indf, mdt);
  g_free(indf);
  assert(indx >= 0 && indx < mdt->nelems);
  *bin_index = indx;
  return TRUE;
}

/** Do some basic setup of an MDT's features. Return TRUE on success. */
gboolean mdt_setup(struct mdt *mdt, const struct mdt_library *mlib,
                   GError **err)
{
  int n, naa, i, ierr;

  /* Do Modeller-specific setup */
  mod_mdt_setup_check(&mdt->base, &mlib->base, &ierr);
  if (ierr != 0) {
    handle_modeller_error(err);
    return FALSE;
  }

  n = naa = 0;
  for (i = 0; i < mdt->base.nfeat; i++) {
    int ifeat = mdt->base.features[i].ifeat;
    const struct mod_mdt_libfeature *feat = &mlib->base.features[ifeat - 1];
    /* you will compare proteins, residues or residue pairs if at least one
       feature requires comparison of proteins, residues or residue pairs,
       respectively. */
    n = MAX(n, feat->iresfeat);
    /* you will generate all NxN residue pairs if any of the relationships
       is asymmetric */
    naa = MAX(naa, feat->isymm);
  }
  mdt->scantype = n;
  mdt->symmetric = (naa == 0);
  return TRUE;
}

/** Get the HDF5 datatype for this MDT. */
hid_t mdt_get_hdf5_type(const struct mod_mdt *mdt)
{
  switch(mdt->bin_type) {
  case MOD_MDTB_FLOAT:
    return H5T_NATIVE_FLOAT;
  case MOD_MDTB_DOUBLE:
    return H5T_NATIVE_DOUBLE;
  case MOD_MDTB_INT32:
    return H5T_NATIVE_INT32;
  case MOD_MDTB_UINT32:
    return H5T_NATIVE_UINT32;
  case MOD_MDTB_INT16:
    return H5T_NATIVE_INT16;
  case MOD_MDTB_UINT16:
    return H5T_NATIVE_UINT16;
  case MOD_MDTB_INT8:
    return H5T_NATIVE_INT8;
  case MOD_MDTB_UINT8:
    return H5T_NATIVE_UINT8;
   default:
     g_assert_not_reached();
     return 0;
   }
}
