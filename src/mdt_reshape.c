/** \file mdt_reshape.c    Functions to reshape MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Set new MDT indices using new offset and shape */
static void reshape_mdt_indices(const struct mod_mdt *mdtin,
                                struct mod_mdt *mdtout, const int offset[],
                                const int shape[], const int features[],
                                const int old_position[])
{
  int i;

  mdtout->nelems = 1;
  for (i = 0; i < mdtin->nfeat; i++) {
    const struct mod_mdt_feature *infeat = &mdtin->features[old_position[i]];
    struct mod_mdt_feature *feat = &mdtout->features[i];
    feat->istart = offset[i] + 1;
    if (shape[i] <= 0) {
      feat->iend = infeat->iend + shape[i];
    } else {
      feat->iend = offset[i] + shape[i];
    }

    feat->ifeat = features[i];
    feat->nbins = feat->iend - feat->istart + 1;
    mdtout->nelems *= feat->nbins;
  }

  make_mdt_stride(mdtout);
}


/** Do the hard work of reshaping the table. */
static void reshape_mdt_table(const struct mod_mdt *mdtin,
                              struct mod_mdt *mdtout, const int new_position[])
{
  int *out_indf, *in_indf;
  out_indf = mdt_start_indices(mdtout);
  in_indf = g_malloc(sizeof(int) * mdtin->nfeat);
  do {
    int i, i1, i2;
    for (i = 0; i < mdtin->nfeat; i++) {
      in_indf[i] = out_indf[new_position[i]];
    }
    i1 = indmdt(in_indf, mdtin);
    i2 = indmdt(out_indf, mdtout);
    mod_mdt_bin_set(mdtout, i2, mod_mdt_bin_get(mdtin, i1));
  } while (roll_ind_mdt(out_indf, mdtout, mdtout->nfeat));
  g_free(in_indf);
  g_free(out_indf);
}

/** Get mapping from old to new features */
static gboolean get_position_mappings(const struct mod_mdt *mdt,
                                      const int features[], int old_position[],
                                      int new_position[], const char *routine,
                                      GError **err)
{
  int i, j;
  for (i = 0; i < mdt->nfeat; i++) {
    int match = 0;
    for (j = 0; j < mdt->nfeat; j++) {
      if (features[i] == mdt->features[j].ifeat) {
        old_position[i] = j;
        new_position[j] = i;
        match = 1;
        break;
      }
    }
    if (!match) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                  "%s: Feature type %d does not exist in input MDT.",
                  routine, features[i]);
      return FALSE;
    }
  }
  return TRUE;
}


/** Check new offset and shape */
static gboolean check_start_end(const struct mod_mdt *mdt, const int offset[],
                                const int shape[], const int old_position[],
                                const int features[], const char *routine,
                                GError **err)
{
  int i;
  for (i = 0; i < mdt->nfeat; i++) {
    const struct mod_mdt_feature *oldfeat = &mdt->features[old_position[i]];
    int end;
    if (shape[i] <= 0) {
      end = oldfeat->iend + shape[i];
    } else {
      end = offset[i] + shape[i];
    }
    if (offset[i] + 1 < oldfeat->istart || end < offset[i] + 1
        || end > oldfeat->iend) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                  "%s: For feature %d, new start %d and size %d are out "
                  "of range.", routine, features[i], offset[i], shape[i]);
      return FALSE;
    }
  }
  return TRUE;
}


/** Reshape an MDT. Return TRUE on success. */
gboolean mdt_reshape(const struct mdt *mdtin, struct mdt *mdtout,
                     const int features[], int n_features, const int offset[],
                     int n_offset, const int shape[], int n_shape,
                     GError **err)
{
  const char *routine = "mdt_reshape";
  int *old_position, *new_position;

  if (n_features != mdtin->base.nfeat || n_offset != mdtin->base.nfeat
      || n_shape != mdtin->base.nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: features, offset and shape must all match"
                " the dimension of the MDT (%d)", routine, mdtin->base.nfeat);
    return FALSE;
  }

  old_position = g_malloc(sizeof(int) * mdtin->base.nfeat);
  new_position = g_malloc(sizeof(int) * mdtin->base.nfeat);

  if (!get_position_mappings(&mdtin->base, features, old_position,
                             new_position, routine, err)
      || !check_start_end(&mdtin->base, offset, shape, old_position, features,
                          routine, err)) {
    g_free(old_position);
    g_free(new_position);
    return FALSE;
  }

  mdt_copy(mdtin, mdtout, mdtin->base.bin_type);
  reshape_mdt_indices(&mdtin->base, &mdtout->base, offset, shape, features,
                      old_position);

  /* reshape the MDT table: */
  reshape_mdt_table(&mdtin->base, &mdtout->base, new_position);

  /* a little heuristic here: */
  if (!mdtout->pdf) {
    mdtout->sample_size = get_mdt_sum(&mdtout->base);
  }

  g_free(old_position);
  g_free(new_position);
  return TRUE;
}
