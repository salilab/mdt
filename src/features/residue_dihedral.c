/** \file residue_dihedral.c  Residue dihedral features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  int dihtype = GPOINTER_TO_INT(data);
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float2_get(&s->dih, residue, dihtype);
  return iclsbin(f, feat);
}

static int add_feature(struct mdt_library *mlib, int protein, int delta,
                       gboolean pos2, const char *name,
                       mdt_dihedral_type dihtype, GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_add(mlib, name, MOD_MDTC_NONE, protein, delta,
                                  pos2, getbin, GINT_TO_POINTER(dihtype), err);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_DIHEDRALS);
  return ifeat;
}

int mdt_feature_chi1_dihedral(struct mdt_library *mlib, int protein,
                              int delta, gboolean pos2, GError **err)
{
  return add_feature(mlib, protein, delta, pos2, "Residue chi1 dihedral",
                     MDT_DIHEDRAL_CHI1, err);
}

int mdt_feature_chi2_dihedral(struct mdt_library *mlib, int protein,
                              int delta, gboolean pos2, GError **err)
{
  return add_feature(mlib, protein, delta, pos2, "Residue chi2 dihedral",
                     MDT_DIHEDRAL_CHI2, err);
}

int mdt_feature_chi3_dihedral(struct mdt_library *mlib, int protein,
                              int delta, gboolean pos2, GError **err)
{
  return add_feature(mlib, protein, delta, pos2, "Residue chi3 dihedral",
                     MDT_DIHEDRAL_CHI3, err);
}
