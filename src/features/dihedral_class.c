/** \file dihedral_class.c  Dihedral class features.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
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
  return mod_int2_get(&s->idihc, residue, dihtype);
}

static int add_feature(struct mdt_library *mlib, int protein, int delta,
                       int align_delta, gboolean pos2, const char *name,
                       mod_mdt_calc precalc, mdt_dihedral_type dihtype,
                       const struct mod_libraries *libs, GError **err)
{
  int ifeat, i;
  struct mod_mdt_libfeature *feat;
  ifeat = mdt_feature_residue_add(mlib, name, precalc, protein, delta,
                                  align_delta, pos2, -1, getbin,
                                  GINT_TO_POINTER(dihtype), NULL, err);
  if (ifeat < 0) {
    return ifeat;
  }
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_DIHEDRALS);

  /* Set number of bins to match number of dihedral classes (plus undef) */
  feat = &mlib->base.features[ifeat - 1];
  mod_mdt_libfeature_nbins_set(feat, libs->resdih.ndihc + 1);
  for (i = 0; i < libs->resdih.ndihc; ++i) {
    g_free(feat->bins[i].symbol);
    feat->bins[i].symbol = g_strdup_printf("%d", i + 1);
    feat->bins[i].rang1 = i;
    feat->bins[i].rang2 = i + 1;
  }
  g_free(feat->bins[libs->resdih.ndihc].symbol);
  feat->bins[libs->resdih.ndihc].symbol = g_strdup("U");
  feat->bins[libs->resdih.ndihc].rang1 = 0;
  feat->bins[libs->resdih.ndihc].rang2 = 0;
  return ifeat;
}

int mdt_feature_chi1_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi1 class", MOD_MDTC_CHI1CL, MDT_DIHEDRAL_CHI1,
                     libs, err);
}

int mdt_feature_chi2_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi2 class", MOD_MDTC_CHI2CL, MDT_DIHEDRAL_CHI2,
                     libs, err);
}

int mdt_feature_chi3_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi3 class", MOD_MDTC_CHI3CL, MDT_DIHEDRAL_CHI3,
                     libs, err);
}

int mdt_feature_chi4_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi4 class", MOD_MDTC_CHI4CL, MDT_DIHEDRAL_CHI4,
                     libs, err);
}

int mdt_feature_chi5_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue chi5 class", MOD_MDTC_CHI5CL, MDT_DIHEDRAL_CHI5,
                     libs, err);
}

int mdt_feature_phi_class(struct mdt_library *mlib, int protein,
                          int delta, int align_delta, gboolean pos2,
                          const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue phi class", MOD_MDTC_PHICL, MDT_DIHEDRAL_PHI,
                     libs, err);
}

int mdt_feature_psi_class(struct mdt_library *mlib, int protein,
                          int delta, int align_delta, gboolean pos2,
                          const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue psi class", MOD_MDTC_PSICL, MDT_DIHEDRAL_PSI,
                     libs, err);
}

int mdt_feature_omega_class(struct mdt_library *mlib, int protein,
                            int delta, int align_delta, gboolean pos2,
                            const struct mod_libraries *libs, GError **err)
{
  return add_feature(mlib, protein, delta, align_delta, pos2,
                     "Residue omega class", MOD_MDTC_OMEGACL,
                     MDT_DIHEDRAL_OMEGA, libs, err);
}
