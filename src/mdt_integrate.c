/** \file mdt_integrate.c  Functions to integrate MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <stdlib.h>
#include <assert.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get the MDT array positions of new features and integrated features */
static gboolean get_feature_indices(const struct mod_mdt *mdt,
                                    const int features[], int n_features,
                                    int **inew_features, int *n_int_features,
                                    int **int_features, const char *routine,
                                    GError **err)
{
  int i, j, indfeat;
  *n_int_features = mdt->nfeat - n_features;
  *inew_features = g_malloc(sizeof(int) * n_features);
  *int_features = g_malloc(sizeof(int) * (*n_int_features));
  for (i = 0; i < n_features; i++) {
    int match = 0;
    for (j = 0; j < mdt->nfeat; j++) {
      if (features[i] == mdt->features[j].ifeat) {
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
      if (mdt->features[i].ifeat == features[j]) {
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
static void copy_mdt_indices_subset(const struct mod_mdt *mdtin,
                                    struct mod_mdt *mdtout,
                                    const int inew_features[])
{
  int i;

  mdtout->nelems = 1;
  for (i = 0; i < mdtout->nfeat; i++) {
    const struct mod_mdt_feature *featin = &mdtin->features[inew_features[i]];
    struct mod_mdt_feature *featout = &mdtout->features[i];

    featout->ifeat = featin->ifeat;
    featout->istart = featin->istart;
    featout->iend = featin->iend;
    featout->nbins = featin->nbins;
    mdtout->nelems *= featout->nbins;
  }
  make_mdt_stride(mdtout);
}

/** Do the actual work of integrating the MDT */
static void integrate_mdt_table(const struct mod_mdt *mdtin,
                                struct mod_mdt *mdtout, int n_features,
                                const int inew_features[], int n_int_features,
                                const int int_features[])
{
  int i, *out_indf, *in_indf;
  in_indf = g_malloc(sizeof(int) * mdtin->nfeat);
  out_indf = mdt_start_indices(mdtout);

  do {
    int i2;
    /* set/reset the indices of the integrated features to istart: */
    for (i = 0; i < n_int_features; i++) {
      in_indf[int_features[i]] = mdtin->features[int_features[i]].istart;
    }
    /* set the indices of the non-integrated features to those in mdtout: */
    for (i = 0; i < n_features; i++) {
      in_indf[inew_features[i]] = out_indf[i];
    }

    /* will be adding the integrated mdtin to this element of mdtout next: */
    i2 = indmdt(out_indf, mdtout);
    mod_mdt_bin_set(mdtout, i2, 0.0);

    /* integrate over all needed dimensions in the first table */
    do {
      int i1 = indmdt(in_indf, mdtin);
      mod_mdt_bin_set(mdtout, i2, mod_mdt_bin_get(mdtout, i2)
                                  + mod_mdt_bin_get(mdtin, i1));

      /* roll the indices of the "integrated" features one forward: */
    } while (roll_inds(in_indf, mdtin, int_features, n_int_features));

    /* roll the indices of the "non-integrated" features one forward: */
  } while (roll_ind_mdt(out_indf, mdtout, mdtout->nfeat));

  g_free(out_indf);
  g_free(in_indf);
}

/** Integrate an MDT. Return TRUE on success. */
gboolean mdt_integrate(const struct mdt *mdtin, struct mdt *mdtout,
                       const int features[], int n_features, GError **err)
{
  static const char *routine = "mdt_integrate";
  int *inew_features, *int_features, n_int_features;

  if (n_features <= 0 || n_features >= mdtin->base.nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: Number of features to be integrated (%d)"
                " must be less than the current number of"
                " features in the MDT (%d) and greater than zero.",
                routine, n_features, mdtin->base.nfeat);
    return FALSE;
  }

  if (!get_feature_indices(&mdtin->base, features, n_features, &inew_features,
                           &n_int_features, &int_features, routine, err)) {
    return FALSE;
  }

  mdt_copy(mdtin, mdtout, mdtin->base.bin_type);
  mdtout->base.nfeat = n_features;
  copy_mdt_indices_subset(&mdtin->base, &mdtout->base, inew_features);

  /* Integrate the MDT table */
  integrate_mdt_table(&mdtin->base, &mdtout->base, n_features, inew_features,
                      n_int_features, int_features);

  /* a little heuristic here */
  if (!mdtout->pdf) {
    mdtout->sample_size = get_mdt_sum(&mdtout->base);
  }
  g_free(inew_features);
  g_free(int_features);
  return TRUE;
}
