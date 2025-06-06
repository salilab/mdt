/** \file residue_dihedral.c  Residue dihedral features.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  int dihtype = GPOINTER_TO_INT(feat->data);
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float2_get(&s->dih, residue, dihtype);
  return feat_to_bin(f, feat);
}

static int add_feature(struct mdt_library *mlib, int protein, int delta,
                       int align_delta, gboolean pos2, const char *name,
                       mdt_dihedral_type dihtype, GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_add(mlib, name, MOD_MDTC_NONE, protein, delta,
                                  align_delta, pos2, -1, getbin,
                                  GINT_TO_POINTER(dihtype), NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_DIHEDRALS);
    mdt_feature_periodic_set(mlib, ifeat, TRUE);
  }
  return ifeat;
}

int mdt_feature_chi1_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi1 dihedral", MDT_DIHEDRAL_CHI1, err);
}

int mdt_feature_chi2_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi2 dihedral", MDT_DIHEDRAL_CHI2, err);
}

int mdt_feature_chi3_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi3 dihedral", MDT_DIHEDRAL_CHI3, err);
}

int mdt_feature_chi4_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi4 dihedral", MDT_DIHEDRAL_CHI4, err);
}

int mdt_feature_phi_dihedral(struct mdt_library *mlib, int protein,
                             int delta, int align_delta, gboolean pos2,
                             GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue phi dihedral", MDT_DIHEDRAL_PHI, err);
}

int mdt_feature_psi_dihedral(struct mdt_library *mlib, int protein,
                             int delta, int align_delta, gboolean pos2,
                             GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue psi dihedral", MDT_DIHEDRAL_PSI, err);
}

int mdt_feature_omega_dihedral(struct mdt_library *mlib, int protein,
                               int delta, int align_delta, gboolean pos2,
                               GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue omega dihedral", MDT_DIHEDRAL_OMEGA, err);
}

int mdt_feature_alpha_dihedral(struct mdt_library *mlib, int protein,
                               int delta, int align_delta, gboolean pos2,
                               GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue alpha dihedral", MDT_DIHEDRAL_ALPHA, err);
}
