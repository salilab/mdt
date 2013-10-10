/** \file xray.c     Protein X-ray resolution feature.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_hdf5.h"

struct feature_data {
  float nmr;
};

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct feature_data *feat_data = (struct feature_data *)feat->data;
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);
  /* artificially change the resolution of the NMR structures
     from the defined -1.00, to decrease the number of
     bins required to hold all defined resolutions while still
     separating NMR from X-ray structures: */
  float f = (seq->resol == -1.00 ? feat_data->nmr : seq->resol);
  return feat_to_bin(f, feat);
}

static gboolean writefunc(hid_t loc_id, const struct mdt_feature *feat,
                          const struct mdt_library *mlib)
{
  struct feature_data *feat_data = (struct feature_data *)feat->data;
  return mdt_hdf5_write_float_attr(loc_id, "nmr", 1, &feat_data->nmr);
}

int mdt_feature_xray_resolution(struct mdt_library *mlib, int protein,
                                float nmr, GError **err)
{
  struct feature_data *feat_data;
  int ifeat;

  feat_data = g_malloc(sizeof(struct feature_data));
  feat_data->nmr = nmr;
  ifeat = mdt_feature_protein_add(mlib, "X-ray resolution", MOD_MDTC_NONE,
                                  protein, getbin, feat_data, g_free, err);
  if (ifeat < 0) {
    g_free(feat_data);
  } else {
    mdt_feature_set_write_callback(mlib, ifeat, writefunc);
  }
  return ifeat;
}
