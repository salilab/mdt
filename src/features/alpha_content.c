/** \file alpha_content.c  Protein alpha content feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  int i, *imnchw, nalpha = 0;
  float f;
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);

  imnchw = mod_int1_pt(&s->imnchw);
  for (i = 0; i < seq->nres; ++i) {
    if (imnchw[i] == 1) {
      nalpha++;
    }
  }
  f = (float)nalpha;
  if (seq->nres > 0) {
    f /= seq->nres;
  }
  return iclsbin(f, feat);
}

int mdt_feature_alpha_content(struct mdt_library *mlib, int protein,
                              GError **err)
{
  int ifeat;
  ifeat = mdt_feature_protein_add(mlib, "Protein alpha content",
                                  MOD_MDTC_MNRAMA, protein, getbin, NULL, err);
  if (ifeat < 0) {
    return ifeat;
  }
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_DIHEDRALS);
  return ifeat;
}
