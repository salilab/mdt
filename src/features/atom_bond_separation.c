/** \file atom_bond_separation.c     Atom-atom bond separation feature.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"
#include "../mdt_residue_bonds.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  int atom1, int atom2, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *struc = mod_alignment_structure_get(aln, protein);
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);
  property_resbond_attyp(aln, protein, prop, mlib, libs, 1);

  int bond_separation = mdt_get_bond_separation(struc, seq, atom1, atom2, prop,
                                                protein,
                                                &mlib->residue_bond_list);
  if (bond_separation == -1) {
    return mdt_feature_undefined_bin_get(feat);
  } else {
    return feat_to_bin(bond_separation, feat);
  }
}

int mdt_feature_atom_bond_separation(struct mdt_library *mlib)
{
  /* Make sure that the residue bonds list is populated */
  mdt_fill_residue_bonds(&mlib->residue_bond_list, mlib, mlib->libs);

  return mdt_feature_atom_pair_add(mlib, "Atom-atom bond separation",
                                   MOD_MDTC_NONE, FALSE, getbin, NULL, NULL);
}
