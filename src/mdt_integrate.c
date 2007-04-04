/** \file mdt_integrate.c  Functions to integrate MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include <assert.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get the MDT array positions of new features and integrated features */
static gboolean get_feature_indices(const struct mdt_type *mdt,
                                    const int features[], int n_features,
                                    int **inew_features, int *n_int_features,
                                    int **int_features, const char *routine,
                                    GError **err)
{
  int i, j, indfeat, *ifeat;
  *n_int_features = mdt->nfeat - n_features;
  *inew_features = g_malloc(sizeof(int) * n_features);
  *int_features = g_malloc(sizeof(int) * (*n_int_features));
  ifeat = f_int1_pt(&mdt->ifeat);
  for (i = 0; i < n_features; i++) {
    int match = 0;
    for (j = 0; j < mdt->nfeat; j++) {
      if (features[i] == ifeat[j]) {
        (*inew_features)[i] = j;
        match = 1;
        break;
      }
    }
    if (!match) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                  "%s: Feature type %d does not exist in input MDT.", routine,
                  features[i]);
      g_free(*inew_features);
      g_free(*int_features);
      return FALSE;
    }
  }

  /* Get the number of MDT array positions of the integrated features
     (complement of new features) */
  indfeat = 0;
  for (i = 0; i < mdt->nfeat; i++) {
    int match = 0;
    for (j = 0; j < n_features; j++) {
      if (ifeat[i] == features[j]) {
        match = 1;
        break;
      }
    }
    if (!match) {
      assert(indfeat < *n_int_features);
      (*int_features)[indfeat++] = i;
    }
  }
  return TRUE;
}

/** Copy a subset of the mdtin feature indices to mdtout */
static void copy_mdt_indices_subset(const struct mdt_type *mdtin,
                                    struct mdt_type *mdtout,
                                    const int inew_features[])
{
  int i;
  int *ifeat_in = f_int1_pt(&mdtin->ifeat);
  int *istart_in = f_int1_pt(&mdtin->istart);
  int *iend_in = f_int1_pt(&mdtin->iend);
  int *nbins_in = f_int1_pt(&mdtin->nbins);
  int *ifeat_out = f_int1_pt(&mdtout->ifeat);
  int *istart_out = f_int1_pt(&mdtout->istart);
  int *iend_out = f_int1_pt(&mdtout->iend);
  int *nbins_out = f_int1_pt(&mdtout->nbins);

  mdtout->nelems = 1;
  for (i = 0; i < mdtout->nfeat; i++) {
    ifeat_out[i] = ifeat_in[inew_features[i]];
    istart_out[i] = istart_in[inew_features[i]];
    iend_out[i] = iend_in[inew_features[i]];
    nbins_out[i] = nbins_in[inew_features[i]];
    mdtout->nelems *= nbins_out[i];
  }
  make_mdt_stride(mdtout);
}

/** Do the actual work of integrating the MDT */
static void integrate_mdt_table(const struct mdt_type *mdtin,
                                struct mdt_type *mdtout, int n_features,
                                const int inew_features[], int n_int_features,
                                const int int_features[])
{
  int i, *out_indf, *in_indf, *in_istart;
  double *out_bin, *in_bin;
  in_indf = g_malloc(sizeof(int) * mdtin->nfeat);
  out_indf = mdt_start_indices(mdtout);
  in_istart = f_int1_pt(&mdtin->istart);
  out_bin = f_double1_pt(&mdtout->bin);
  in_bin = f_double1_pt(&mdtin->bin);

  do {
    int i2;
    /* set/reset the indices of the integrated features to istart: */
    for (i = 0; i < n_int_features; i++) {
      in_indf[int_features[i]] = in_istart[int_features[i]];
    }
    /* set the indices of the non-integrated features to those in mdtout: */
    for (i = 0; i < n_features; i++) {
      in_indf[inew_features[i]] = out_indf[i];
    }

    /* will be adding the integrated mdtin to this element of mdtout next: */
    i2 = indmdt(out_indf, mdtout);
    out_bin[i2] = 0.0;

    /* integrate over all needed dimensions in the first table */
    do {
      int i1 = indmdt(in_indf, mdtin);
      out_bin[i2] += in_bin[i1];

    /* roll the indices of the "integrated" features one forward: */
    } while (roll_inds(in_indf, f_int1_pt(&mdtin->istart),
                       f_int1_pt(&mdtin->iend), mdtin->nfeat,
                       int_features, n_int_features));

  /* roll the indices of the "non-integrated" features one forward: */
  } while (roll_ind(out_indf, f_int1_pt(&mdtout->istart),
                    f_int1_pt(&mdtout->iend), mdtout->nfeat));

  g_free(out_indf);
  g_free(in_indf);
}

/** Integrate an MDT. Return TRUE on success. */
gboolean mdt_integrate(const struct mdt_type *mdtin, struct mdt_type *mdtout,
                       const int features[], int n_features, GError **err)
{
  static const char *routine = "mdt_integrate";
  int *inew_features, *int_features, n_int_features;

  if (n_features <= 0 || n_features >= mdtin->nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: Number of features to be integrated (%d)"
                " must be less than the current number of"
                " features in the MDT (%d) and greater than zero.",
                routine, n_features, mdtin->nfeat);
    return FALSE;
  }

  if (!get_feature_indices(mdtin, features, n_features, &inew_features,
                           &n_int_features, &int_features, routine, err)) {
    return FALSE;
  }

  copy_mdt(mdtin, mdtout);
  mdtout->nfeat = n_features;
  copy_mdt_indices_subset(mdtin, mdtout, inew_features);

  /* Integrate the MDT table */
  integrate_mdt_table(mdtin, mdtout, n_features, inew_features, n_int_features,
                      int_features);

  /* a little heuristic here */
  if (!mdtout->pdf) {
    mdtout->sample_size = get_sum(f_double1_pt(&mdtout->bin), mdtout->nelems);
  }
  g_free(inew_features);
  g_free(int_features);
  return TRUE;
}
