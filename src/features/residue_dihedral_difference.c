/** \file residue_dihedral_difference.c  Residue dihedral difference features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"

static int getbin(const struct mod_alignment *aln, int protein1, int protein2,
                  int alnpos, struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  int dihtype = GPOINTER_TO_INT(data);
  struct mod_structure *s1 = mod_alignment_structure_get(aln, protein1);
  struct mod_structure *s2 = mod_alignment_structure_get(aln, protein2);
  int residue1 = mod_int2_get(&aln->ialn, alnpos, protein1) - 1;
  int residue2 = mod_int2_get(&aln->ialn, alnpos, protein2) - 1;
  float d1 = mod_float2_get(&s1->dih, residue1, dihtype);
  float d2 = mod_float2_get(&s2->dih, residue2, dihtype);
  float diff;

  /* If either dihedral is out of range (-999), return undefined */
  if (d1 == -999 || d2 == -999) {
    return feat->nbins;
  }

  /* Take difference and ensure it's in the range -180 to 180 */
  diff = d2 - d1;
  if (diff < -180.0) {
    diff += 360.0;
  }
  if (diff > 180.0) {
    diff -= 360.0;
  }
  return iclsbin(diff, feat);
}

static int add_feature(struct mdt_library *mlib, int protein1, int protein2,
                       const char *name, mdt_dihedral_type dihtype,
                       GError **err)
{
  int ifeat;
  ifeat = mdt_feature_aligned_residue_add(mlib, name, MOD_MDTC_NONE, protein1,
                                          protein2, getbin,
                                          GINT_TO_POINTER(dihtype), err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_DIHEDRALS);
  }
  return ifeat;
}

int mdt_feature_phi_dihedral_difference(struct mdt_library *mlib, int protein1,
                                        int protein2, GError **err)
{
  return add_feature(mlib, protein1, protein2, "Phi dihedral difference",
                     MDT_DIHEDRAL_PHI, err);
}

int mdt_feature_psi_dihedral_difference(struct mdt_library *mlib, int protein1,
                                        int protein2, GError **err)
{
  return add_feature(mlib, protein1, protein2, "Psi dihedral difference",
                     MDT_DIHEDRAL_PSI, err);
}

int mdt_feature_omega_dihedral_difference(struct mdt_library *mlib,
                                          int protein1, int protein2,
                                          GError **err)
{
  return add_feature(mlib, protein1, protein2, "Omega dihedral difference",
                     MDT_DIHEDRAL_OMEGA, err);
}