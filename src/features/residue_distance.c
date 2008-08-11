/** \file residue_distance.c     Residue-residue distance feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue1,
                  int residue2, struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  int atom1 = mod_int1_get(&s->idsta1, residue1) - 1;
  int atom2 = mod_int1_get(&s->idsta1, residue2) - 1;
  return idist0(atom1, atom2, s, feat);
}

int mdt_feature_residue_distance(struct mdt_library *mlib, int protein,
                                 GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_pair_add(mlib, "Residue-residue distance",
                                       MOD_MDTC_ATIND, protein, TRUE, getbin,
                                       NULL, err);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}
