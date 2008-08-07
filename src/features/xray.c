/** \file xray.c     Protein X-ray resolution feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);
  /* artificially change the resolution of the NMR structures
     from the defined -1.00 to 0.45, to decrease the number of
     bins required to hold all defined resolutions while still
     separating NMR from X-ray structures: */
  float f = (seq->resol == -1.00 ? 0.45 : seq->resol);
  return iclsbin(f, feat);
}

int mdt_feature_xray_resolution(struct mdt_library *mlib, int protein,
                                GError **err)
{
  return mdt_feature_protein_add(mlib, "X-ray resolution", MOD_MDTC_NONE,
                                 protein, getbin, NULL, err);
}
