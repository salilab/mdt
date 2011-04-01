/** \file mdt_super_smooth.c  Functions to smooth MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "modeller.h"
#include "mdt.h"
#include "num_recipes.h"
#include "util.h"

/** All information needed for a single combination of features */
struct combination {
  /** MDT stride for previous level */
  int *stride1;
  /** MDT stride for current level */
  int *stride2;
  /** Number of bins in MDT at this level */
  int *ndims2;
  /** Indices of fixed features for the previous level */
  int *i_feat_fixn1;
  /** Indices of fixed features for the current level */
  int *i_feat_fixn2;
  /** Start index into total bin array for this level */
  int i1add;
  /** End index into total bin array for this level */
  int i2add;
};

/** A vector of combination structures */
struct combination_vector {
  /** The combination data */
  struct combination *combinations;
  /** Number of combination structs */
  int length;
  /** Allocated number of combination structs */
  int allocated;
};

/** Create and return a new combination_vector */
static struct combination_vector *new_combination_vector(void)
{
  struct combination_vector *vec = g_malloc(sizeof(struct combination_vector));
  vec->combinations = NULL;
  vec->length = vec->allocated = 0;
  return vec;
}

/** Free memory used by a combination vector */
static void free_combination_vector(struct combination_vector *vec)
{
  int i;
  for (i = 0; i < vec->allocated; i++) {
    g_free(vec->combinations[i].stride1);
    g_free(vec->combinations[i].stride2);
    g_free(vec->combinations[i].ndims2);
    g_free(vec->combinations[i].i_feat_fixn1);
    g_free(vec->combinations[i].i_feat_fixn2);
  }
  g_free(vec->combinations);
  g_free(vec);
}

/** Calculate number of feature combinations, and reallocate arrays if
    necessary. */
static void calculate_combinations(struct combination_vector *vec,
                                   const struct mod_mdt *mdt, int dimensions,
                                   int n_feat_fix)
{
  vec->length = nperm(mdt->nfeat - dimensions)
      / (nperm(n_feat_fix) * nperm(mdt->nfeat - dimensions - n_feat_fix));
  if (vec->length > vec->allocated) {
    int i;
    int maxcom = MAX(vec->allocated * 3 / 2, vec->length);
    maxcom = MAX(maxcom, 10);
    vec->combinations = g_realloc(vec->combinations,
                                  sizeof(struct combination) * maxcom);
    for (i = vec->allocated; i < maxcom; i++) {
      struct combination *com = &vec->combinations[i];
      com->stride1 = g_malloc(sizeof(int) * mdt->nfeat);
      com->stride2 = g_malloc(sizeof(int) * mdt->nfeat);
      com->ndims2 = g_malloc(sizeof(int) * mdt->nfeat);
      com->i_feat_fixn1 = g_malloc(sizeof(int) * mdt->nfeat);
      com->i_feat_fixn2 = g_malloc(sizeof(int) * mdt->nfeat);
    }
    vec->allocated = maxcom;
  }
}

/** Normalize an array of frequencies, and return its sum */
static double normalize_freq(double frq[], int nbin)
{
  static const float divisor = 1e-15;
  double sumfrq = get_sum(frq, nbin);
  if (sumfrq > divisor) {
    int i;
    for (i = 0; i < nbin; i++) {
      frq[i] /= sumfrq;
    }
  } else {
    int i;
    for (i = 0; i < nbin; i++) {
      frq[i] = 1.0 / (double)nbin;
    }
  }
  return sumfrq;
}


/** Returns p^n of Eq. 10 */
static void smthfrq(const double apriori[], double frq[], int i_val_fix[],
                    const int stride2[], const int i_start_fix[], int level,
                    double prior_weight, int nbinx, double bin2[], int nelems)
{
  int i, is1;
  float sumfrq, w1, w2;

  /* sum(apriori(1:nbinx)) can be zero when no smoothing is done
     (ie, prior_weight = 0). w1 will not be zero, and that is
     ok because w1*apriori is still 0. */

  /* w2 may also be 0 if there are no data points in frq, but then frq
     will be normalized to the even distribution by normalize_freq,
     so all is ok, even when apriori is 0. */

  /* sumfrq can be zero when there are no points in the data, in which
     case w2 = 0 and pdf frq(1:nbinx) does not matter: */

  sumfrq = normalize_freq(frq, nbinx);

  /* the weights for the 'a priori' and 'data' pdf's (w1 + w2 = 1) */
  weights(prior_weight, nbinx, sumfrq, &w1, &w2);

  /* smooth: */
  i_val_fix[level - 1] = 1;
  is1 = indmdt_full(i_val_fix, stride2, level, i_start_fix);
  for (i = 0; i < nbinx; i++) {
    bin2[is1 + i] = w1 * apriori[i] + w2 * frq[i];
  }

  /* make sure bin2 is nicely normalized: */
  normalize_freq(&bin2[is1], nbinx);
}


/** If combination 'comb' is a subset of combination ic2, update inds1 and
    return true; otherwise, return false. */
static gboolean get_inds_combination(const struct combination *comb,
                                     int n_feat_fix, const int i_feat_fix[],
                                     const int i_val_fix[], int inds1[])
{
  int i, j;
  for (i = 0; i < n_feat_fix - 1; i++) {
    gboolean success = FALSE;
    for (j = 0; j < n_feat_fix && !success; j++) {
      if (comb->i_feat_fixn1[i] == i_feat_fix[j]) {
        /* successful: feature i from combination i1 occurs in combination
           ic2: remember the index for the feature value */
        inds1[i] = i_val_fix[j];
        success = TRUE;
      }
    }
    if (!success) {
      /* not successful: i1 combination does not occur within ic2
         combination: */
      return FALSE;
    }
  }
  return TRUE;
}


/** Return the a priori distribution, A^n in Eq. 12. */
static void getapriori(gboolean entropy_weighing, const struct mod_mdt *mdt,
                       const struct combination_vector *vec, int nbinx,
                       const int i_feat_fix[], const int i_val_fix[],
                       int ncomb1, int n_feat_fix, const int i_start_fix[],
                       int numb_features, double apriori[])
{
  if (ncomb1 == 0) {
    int i;
    for (i = 0; i < nbinx; i++) {
      apriori[i] = 1.0 / (double)nbinx;
    }
  } else {
    float emax;
    int i, i1, *inds1;

    inds1 = g_malloc(sizeof(int) * n_feat_fix);

    for (i = 0; i < nbinx; i++) {
      apriori[i] = 0.0;
    }
    emax = log((float)nbinx);

    /* for ALL combinations from the previous level: */
    for (i1 = 0; i1 < ncomb1; i1++) {
      const struct combination *comb = &vec->combinations[i1];
      if (get_inds_combination(comb, n_feat_fix, i_feat_fix, i_val_fix,
                               inds1)) {
        float w;
        int is1;
        /* calculate the weight (prop. to negentropy) of the ncmbtst-th
           conditional distribution of x given features in combinations
           1 .. ncomb1; be careful: some pdf's may be 0
           (their weight must be 0 and the other weights larger - better
           solved by normalization at the end)! */

        /* from the previous level: */
        inds1[n_feat_fix - 1] = 1;
        is1 = comb->i1add + indmdt_full(inds1, comb->stride1, n_feat_fix,
                                        i_start_fix);

        if (entropy_weighing) {
          /* w is proportional to rho_c in Eq. 12; w is the nominator in
             Eq. 14, also defined by Eq. 13: */
          double *tmparr = (double *)g_malloc(sizeof(double) * nbinx);
          for (i = 0; i < nbinx; ++i) {
            tmparr[i] = mod_mdt_bin_get(mdt, is1 + i);
          }
          w = emax - entrp1(tmparr, nbinx);
          g_free(tmparr);
        } else {
          w = 1.0;
        }

        for (i = 0; i < nbinx; i++) {
          apriori[i] += w * mod_mdt_bin_get(mdt, is1 + i);
        }
      }
    }

    /* normalize (because weights w do not generally sum to 1): */
    normalize_freq(apriori, nbinx);

    g_free(inds1);
  }
}


/** Finished with level level; copy the new variables to the old ones: */
static void finish_level(int level, int n_feat_fix, int nelm2,
                         struct combination_vector *vec, struct mod_mdt *mdt,
                         int *ncomb1, double bin2[])
{
  int icomb, i;
  *ncomb1 = vec->length;

  for (icomb = 0; icomb < vec->length; icomb++) {
    int ifeat;
    struct combination *comb = &vec->combinations[icomb];

    comb->i1add = comb->i2add;
    for (ifeat = 0; ifeat < level; ifeat++) {
      comb->stride1[ifeat] = comb->stride2[ifeat];
    }
    for (ifeat = 0; ifeat < n_feat_fix; ifeat++) {
      comb->i_feat_fixn1[ifeat] = comb->i_feat_fixn2[ifeat];
    }
  }
  for (i = 0; i < nelm2; ++i) {
    mod_mdt_bin_set(mdt, i, bin2[i]);
  }
}


/** Prepare arrays for smoothing at level level */
static void prepare_level(int level, struct combination_vector *vec,
                          int nbinx, int dimensions, int *nelm2, int *maxelm2,
                          double **bin2, const struct mod_mdt *mdtin,
                          int n_feat_fix)
{
  int *i_feat_fix = NULL, ic2;
  /* get the number of different combinations of n_feat_fix features
     from the set of nfeat-dimensions possible features */
  calculate_combinations(vec, mdtin, dimensions, n_feat_fix);

  /* number of elements in the bin2 array for all 1 .. ncomb2 combinations
     of independent features; bin2 stores ncomb2 MDT's linearly ... */
  *nelm2 = 0;

  /* prepare the indices for counting feature indices at this level
     (save them for the next cycle as well) */
  ic2 = 0;
  while (roll_ind_comb(&i_feat_fix, n_feat_fix, mdtin->nfeat - dimensions)) {
    int i;
    struct combination *comb = &vec->combinations[ic2];

    /* the BIN2 index of the first element of the ic2-th combination: */
    comb->i2add = *nelm2;

    for (i = 0; i < n_feat_fix; i++) {
      comb->ndims2[i] = mdtin->features[i_feat_fix[i]].nbins;
      /* save the feature combination for the next level: */
      comb->i_feat_fixn2[i] = i_feat_fix[i];
    }
    comb->ndims2[level - 1] = nbinx;

    *nelm2 += make_mdt_stride_full(comb->ndims2, level, comb->stride2);
    ic2++;
  }
  g_free(i_feat_fix);

  if (*nelm2 > *maxelm2) {
    *maxelm2 = MAX(*maxelm2 * 3 / 2, *nelm2);
    *bin2 = g_realloc(*bin2, sizeof(double) * (*maxelm2));
  }
}

/** Scan all generated combinations of independent features at level level */
static void build_level_combination(int level, int dimensions, int nbinx,
                                    struct combination_vector *vec,
                                    int icomb, int n_feat_fix,
                                    int i_feat_fix[], int n_bins_fix[],
                                    int i_val_fix[], int i_start_fix[],
                                    const struct mod_mdt *mdtin,
                                    const struct mod_mdt *mdtout,
                                    gboolean entropy_weighing, int ncomb1,
                                    int nelm2, double apriori[], double frq[],
                                    double prior_weight, double bin2[])
{
  int i;
  struct combination *comb = &vec->combinations[icomb];

  for (i = 0; i < n_feat_fix; i++) {
    i_feat_fix[i] = comb->i_feat_fixn2[i];
    n_bins_fix[i] = mdtin->features[i_feat_fix[i]].nbins;
  }
  for (i = 0; i < mdtin->nfeat; i++) {
    i_val_fix[i] = 1;
    i_start_fix[i] = 1;
  }

  /* calculate the new smoothed pdf ic2 at level LEVEL by summing pdfs
     corresponding to combinations of n_feat_fix independent features
     included in level LEVEL */

  /* for each generated combination of features i_feat_fix, scan all
     combinations of values i_val_fix(1:n_feat_fix) of features
     i_feat_fix(1:n_feat_fix): */

  do {
    /* calculate the a priori pdf ic2 for the current i_val_fix,
       apriori[]      A^n         (a priori distribution)
       mdtin->bin[]   W'          (raw frequencies)
       mdtout->bin[]  p^(n-1)     (smoothed distribution at level-1)
       frq[]          W^n         (normalized frequencies) */

    /* initialize APRIORI, based on mdtout->bin (which depends on what happened
       in the previous LEVEL cycle);
       apriori will be 1/nbinx, if prior_weight = 0 and no data:
       this routine does the job of Eq. 12: */
    getapriori(entropy_weighing, mdtout, vec, nbinx,
               i_feat_fix, i_val_fix, ncomb1, n_feat_fix, i_start_fix,
               mdtin->nfeat, apriori);

    /* calc the experimental frequency ic2 for the current i_val_fix:
       have to sum over all values for all independent features, except
       for the n_feat_fix features which are fixed to their current
       values in i_val_fix */

    /* initializes FRQ based on BIN;
       this routine does the job of Eq. 11 (without normalization): */
    getfrq(mdtin, i_feat_fix, n_feat_fix, i_val_fix, dimensions, nbinx, frq);

    /* calculate bin2 = weighted sum of the a priori and experimental
       pdf's: sets BIN2, based on APRIORI and FRQ
       this routine does normalization of Eq. 11 and smoothing of Eq. 10: */
    smthfrq(apriori, frq, i_val_fix, comb->stride2, i_start_fix, level,
            prior_weight, nbinx, &bin2[comb->i2add], nelm2 - comb->i2add);

    /* more combinations of values of independent variables? */
  } while (roll_ind(i_val_fix, i_start_fix, n_bins_fix, n_feat_fix));
}


/** Do one level of smoothing */
static void super_smooth_level(const struct mod_mdt *mdtin,
                               struct mod_mdt *mdtout, float prior_weight,
                               gboolean entropy_weighing, int level, int nbinx,
                               int *ncomb1, struct combination_vector *vec,
                               int dimensions, int i_feat_fix[],
                               int n_bins_fix[], int i_val_fix[],
                               int i_start_fix[], double apriori[],
                               double frq[], double **bin2,
                               int *nelm2, int *maxelm2)
{
  int n_feat_fix, i;

  /* number of fixed features in each combination of independent features
     at level level */
  n_feat_fix = level - 1;

  prepare_level(level, vec, nbinx, dimensions, nelm2, maxelm2, bin2, mdtin,
                n_feat_fix);

  for (i = 0; i < vec->length; i++) {
    build_level_combination(level, dimensions, nbinx, vec, i, n_feat_fix,
                            i_feat_fix, n_bins_fix, i_val_fix, i_start_fix,
                            mdtin, mdtout, entropy_weighing, *ncomb1, *nelm2,
                            apriori, frq, prior_weight, *bin2);
  }

  finish_level(level, n_feat_fix, *nelm2, vec, mdtout, ncomb1, *bin2);
}


/** Get the number of bins for dependent features. Return TRUE on success. */
static gboolean get_nbinx(const struct mdt *mdt, int dimensions, int *nbinx,
                          GError **err)
{
  int nfeat = mdt->base.nfeat;
  if (dimensions <= 0 || dimensions >= nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "'dimensions' is %d; it must be greater than 0, and less than "
                "the dimensionality of this MDT (%d)", dimensions, nfeat);
    return FALSE;
  } else {
    int i;
    *nbinx = 1;
    for (i = 0; i < dimensions; i++) {
      (*nbinx) *= mdt->base.features[nfeat - 1 - i].nbins;
    }
    return TRUE;
  }
}


/** Super-duper multi-level hierarchical recursive multi-dimensional
    smoothing of sparse MDT frequency tables. Return TRUE on success. */
gboolean mdt_super_smooth(const struct mdt *mdtin, struct mdt *mdtout,
                          int dimensions, float prior_weight,
                          gboolean entropy_weighing, GError **err)
{
  int level, ncomb1, *i_feat_fix, *n_bins_fix, *i_val_fix, *i_start_fix,
      nelm2 = 0, maxelm2 = 0, nbinx;
  double *apriori, *frq, *bin2;
  struct combination_vector *vec;

  if (!get_nbinx(mdtin, dimensions, &nbinx, err)) {
    return FALSE;
  }

  vec = new_combination_vector();
  mdt_copy(mdtin, mdtout, mdtin->base.bin_type);

  i_feat_fix = g_malloc(sizeof(int) * mdtin->base.nfeat);
  n_bins_fix = g_malloc(sizeof(int) * mdtin->base.nfeat);
  i_start_fix = g_malloc(sizeof(int) * mdtin->base.nfeat);
  i_val_fix = g_malloc(sizeof(int) * mdtin->base.nfeat);
  apriori = g_malloc(sizeof(double) * nbinx);
  frq = g_malloc(sizeof(double) * nbinx);
  bin2 = NULL;

  /* initialize variables 'from the previous level' */
  ncomb1 = 0;

  for (level = 1; level <= mdtin->base.nfeat - dimensions + 1; level++) {
    super_smooth_level(&mdtin->base, &mdtout->base, prior_weight,
                       entropy_weighing, level, nbinx, &ncomb1, vec, dimensions,
                       i_feat_fix, n_bins_fix, i_val_fix, i_start_fix,
                       apriori, frq, &bin2, &nelm2, &maxelm2);
  }

  free_combination_vector(vec);
  g_free(i_feat_fix);
  g_free(n_bins_fix);
  g_free(i_start_fix);
  g_free(i_val_fix);
  g_free(apriori);
  g_free(frq);
  g_free(bin2);
  return TRUE;
}
