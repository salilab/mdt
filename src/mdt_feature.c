/** \file mdt_feature.c    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "mdt_feature.h"
#include "mdt_atom_classes.h"
#include "mdt_index.h"

void mdt_feature_add_needed_file(struct mdt_library *mlib, int ifeat,
                                 mod_mdt_file filetype)
{
  struct mod_mdt_libfeature *feat = &mlib->base.features[ifeat - 1];
  feat->ndatfil++;
  feat->idatfil = g_realloc(feat->idatfil, sizeof(int) * feat->ndatfil);
  feat->idatfil[feat->ndatfil - 1] = filetype;
}

void mdt_feature_set_write_callback(struct mdt_library *mlib, int ifeat,
                                    mdt_cb_feature_write writefunc)
{
  struct mdt_feature *feat;
  feat = &g_array_index(mlib->features, struct mdt_feature, ifeat - 1);
  feat->writefunc = writefunc;
}

void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins)
{
  struct mod_mdt_libfeature *libfeat;
  libfeat = &mlib->base.features[ifeat - 1];

  /* Add one extra for the undefined bin */
  mod_mdt_libfeature_nbins_set(libfeat, nbins + 1);
  libfeat->bins[nbins].rang1 = 0.;
  libfeat->bins[nbins].rang2 = 1.;
  g_free(libfeat->bins[nbins].symbol);
  libfeat->bins[nbins].symbol = g_strdup("U");
}

void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol)
{
  struct mod_mdt_libfeature *libfeat;
  libfeat = &mlib->base.features[ifeat - 1];

  libfeat->bins[bin].rang1 = start;
  libfeat->bins[bin].rang2 = end;
  g_free(libfeat->bins[bin].symbol);
  libfeat->bins[bin].symbol = g_strdup(symbol ? symbol : "");
}

void mdt_feature_periodic_set(struct mdt_library *mlib, int ifeat,
                              gboolean periodic)
{
  struct mdt_feature *feat;
  feat = &g_array_index(mlib->features, struct mdt_feature, ifeat - 1);
  feat->periodic = periodic;
}

gboolean mdt_feature_periodic_get(const struct mdt_library *mlib, int ifeat)
{
  struct mdt_feature *feat;
  feat = &g_array_index(mlib->features, struct mdt_feature, ifeat - 1);
  return feat->periodic;
}

/** \return TRUE iff 'protein' is in range. */
static gboolean check_protein(int protein, const char *feattype,
                              mod_mdt_protein *mdt_protein, GError **err)
{
  switch(protein) {
  case 0:
    *mdt_protein = MOD_MDTP_A;
    break;
  case 1:
    *mdt_protein = MOD_MDTP_B;
    break;
  case 2:
    *mdt_protein = MOD_MDTP_C;
    break;
  default:
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s features can act only on protein 0, 1, or 2; "
                "%d was given", feattype, protein);
    return FALSE;
  }
  return TRUE;
}

/** \return TRUE iff the protein pair is valid. */
static gboolean check_protein_pair(int protein1, int protein2,
                                   const char *feattype, GError **err)
{
  if (protein1 != 0 || (protein2 != 1 && protein2 != 2)) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s features can act only on protein pairs (0,1) or (0,2);"
                "(%d,%d) was given", feattype, protein1, protein2);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Helper function to add a new feature structure, and return it. */
static struct mdt_feature *add_feature(struct mdt_library *mlib, int *nfeat,
                                       mdt_feature_type type, void *data,
                                       mdt_cb_free freefunc)
{
  struct mdt_feature *newfeat;
  mlib->feature_added = TRUE;
  *nfeat = mlib->base.nfeat + 1;
  mlib->features = g_array_set_size(mlib->features, *nfeat);
  newfeat = &g_array_index(mlib->features, struct mdt_feature, *nfeat - 1);
  newfeat->periodic = FALSE;
  newfeat->type = type;
  newfeat->data = data;
  newfeat->freefunc = freefunc;
  newfeat->writefunc = NULL;
  return newfeat;
}

int mdt_feature_protein_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein,
                            mdt_cb_feature_protein getbin, void *data,
                            mdt_cb_free freefunc, GError **err)
{
  GString *fullname;
  struct mdt_feature_protein *feat;
  mod_mdt_protein mdt_protein;
  int nfeat;

  if (!check_protein(protein, "Protein", &mdt_protein, err)) {
    return -1;
  }

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_PROTEIN, data,
                       freefunc)->u.protein);
  feat->protein = protein;
  feat->getbin = getbin;
  fullname = g_string_new(name);
  g_string_append_printf(fullname, " of protein %d", protein);
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              mdt_protein, MOD_MDTS_PROTEIN, FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}

int mdt_feature_protein_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein1,
                                 int protein2,
                                 mdt_cb_feature_protein_pair getbin, void *data,
                                 mdt_cb_free freefunc, GError **err)
{
  char *fullname;
  struct mdt_feature_protein_pair *feat;
  int nfeat;

  if (!check_protein_pair(protein1, protein2, "Protein pair", err)) {
    return -1;
  }

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_PROTEIN_PAIR,
                       data, freefunc)->u.protein_pair);
  feat->protein1 = protein1;
  feat->protein2 = protein2;
  feat->getbin = getbin;
  fullname = g_strdup_printf("%s of proteins (%d,%d)", name, protein1,
                             protein2);
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname, precalc_type,
                              protein2 == 1 ? MOD_MDTP_AB : MOD_MDTP_AC,
                              MOD_MDTS_PROTEIN, FALSE, 0);
  g_free(fullname);
  return nfeat;
}

int mdt_feature_residue_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein, int delta,
                            int align_delta, gboolean pos2,
                            int bin_seq_outrange, mdt_cb_feature_residue getbin,
                            void *data, mdt_cb_free freefunc, GError **err)
{
  GString *fullname;
  struct mdt_feature_residue *feat;
  mod_mdt_protein mdt_protein;
  int nfeat;

  if (!check_protein(protein, "Residue", &mdt_protein, err)) {
    return -1;
  }

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_RESIDUE, data,
                       freefunc)->u.residue);
  feat->protein = protein;
  feat->delta = delta;
  feat->align_delta = align_delta;
  feat->pos2 = pos2;
  feat->bin_seq_outrange = bin_seq_outrange;
  feat->getbin = getbin;
  fullname = g_string_new(name);
  g_string_append_printf(fullname, " of protein %d", protein);
  if (pos2) {
    g_string_append(fullname, " at pos2");
  }
  if (delta != 0) {
    g_string_append_printf(fullname, ", at delta %d", delta);
  }
  if (align_delta != 0) {
    g_string_append_printf(fullname, ", at alignment delta %d", align_delta);
  }
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              mdt_protein,
                              pos2 ? MOD_MDTS_RESIDUE_PAIR : MOD_MDTS_RESIDUE,
                              FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}

int mdt_feature_residue_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein,
                                 gboolean asymmetric,
                                 mdt_cb_feature_residue_pair getbin, void *data,
                                 mdt_cb_free freefunc, GError **err)
{
  struct mdt_feature_residue_pair *feat;
  mod_mdt_protein mdt_protein;
  int nfeat;

  if (!check_protein(protein, "Residue pair", &mdt_protein, err)) {
    return -1;
  }

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_RESIDUE_PAIR,
                       data, freefunc)->u.residue_pair);
  feat->protein = protein;
  feat->getbin = getbin;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              mdt_protein, MOD_MDTS_RESIDUE_PAIR, asymmetric,
                              0);
  return nfeat;
}

int mdt_feature_aligned_residue_add(struct mdt_library *mlib, const char *name,
                                    mod_mdt_calc precalc_type, int protein1,
                                    int protein2,
                                    mdt_cb_feature_aligned_residue getbin,
                                    void *data, mdt_cb_free freefunc,
                                    GError **err)
{
  char *fullname;
  struct mdt_feature_aligned_residue *feat;
  int nfeat;

  if (!check_protein_pair(protein1, protein2, "Aligned residue", err)) {
    return -1;
  }

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_ALIGNED_RESIDUE,
                       data, freefunc)->u.aligned_residue);
  feat->protein1 = protein1;
  feat->protein2 = protein2;
  feat->getbin = getbin;
  fullname = g_strdup_printf("%s of proteins (%d,%d)", name, protein1,
                             protein2);
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname, precalc_type,
                              protein2 == 1 ? MOD_MDTP_AB : MOD_MDTP_AC,
                              MOD_MDTS_RESIDUE, FALSE, 0);
  g_free(fullname);
  return nfeat;
}

int mdt_feature_aligned_residue_pair_add(
    struct mdt_library *mlib, const char *name, mod_mdt_calc precalc_type,
    int protein1, int protein2, gboolean asymmetric,
    mdt_cb_feature_aligned_residue_pair getbin, void *data,
    mdt_cb_free freefunc, GError **err)
{
  char *fullname;
  struct mdt_feature_aligned_residue_pair *feat;
  int nfeat;

  if (!check_protein_pair(protein1, protein2, "Aligned residue pair", err)) {
    return -1;
  }

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_ALIGNED_RESIDUE_PAIR,
                       data, freefunc)->u.aligned_residue_pair);
  feat->protein1 = protein1;
  feat->protein2 = protein2;
  feat->getbin = getbin;
  fullname = g_strdup_printf("%s of proteins (%d,%d)", name, protein1,
                             protein2);
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname, precalc_type,
                              protein2 == 1 ? MOD_MDTP_AB : MOD_MDTP_AC,
                              MOD_MDTS_RESIDUE_PAIR, asymmetric, 0);
  g_free(fullname);
  return nfeat;
}

int mdt_feature_atom_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type, gboolean pos2,
                         mdt_cb_feature_atom getbin, void *data,
                         mdt_cb_free freefunc)
{
  GString *fullname;
  struct mdt_feature_atom *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_ATOM, data, freefunc)->u.atom);
  feat->pos2 = pos2;
  feat->getbin = getbin;
  fullname = g_string_new(name);
  if (pos2) {
    g_string_append(fullname, " at pos2");
  }
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              MOD_MDTP_A,
                              pos2 ? MOD_MDTS_ATOM_PAIR : MOD_MDTS_ATOM,
                              FALSE, 0);
  g_string_free(fullname, TRUE);
  mdt_feature_add_needed_file(mlib, nfeat, MOD_MDTF_STRUCTURE);
  return nfeat;
}

int mdt_feature_atom_pair_add(struct mdt_library *mlib, const char *name,
                              mod_mdt_calc precalc_type, gboolean asymmetric,
                              mdt_cb_feature_atom_pair getbin, void *data,
                              mdt_cb_free freefunc)
{
  struct mdt_feature_atom_pair *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_ATOM_PAIR, data,
                       freefunc)->u.atom_pair);
  feat->getbin = getbin;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_ATOM_PAIR, asymmetric, 0);
  mdt_feature_add_needed_file(mlib, nfeat, MOD_MDTF_STRUCTURE);
  return nfeat;
}

int mdt_feature_tuple_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, gboolean pos2,
                          mdt_cb_feature_tuple getbin, void *data,
                          mdt_cb_free freefunc)
{
  GString *fullname;
  struct mdt_feature_tuple *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_TUPLE, data,
                       freefunc)->u.tuple);
  feat->pos2 = pos2;
  feat->getbin = getbin;
  fullname = g_string_new(name);
  if (pos2) {
    g_string_append(fullname, " at pos2");
  }
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              MOD_MDTP_A,
                              pos2 ? MOD_MDTS_TUPLE_PAIR : MOD_MDTS_TUPLE,
                              FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}

int mdt_feature_tuple_pair_add(struct mdt_library *mlib, const char *name,
                               mod_mdt_calc precalc_type,
                               mdt_cb_feature_tuple_pair getbin, void *data,
                               mdt_cb_free freefunc)
{
  struct mdt_feature_tuple_pair *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_TUPLE_PAIR,
                       data, freefunc)->u.tuple_pair);
  feat->getbin = getbin;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR, TRUE, 0);
  return nfeat;
}

int mdt_feature_bond_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type,
                         mdt_cb_feature_bond getbin, void *data,
                         mdt_cb_free freefunc)
{
  struct mdt_feature_bond *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_BOND, data, freefunc)->u.bond);
  feat->type = MDT_BOND_TYPE_BOND;
  feat->getbin = getbin;

  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_BOND, FALSE, 0);
  return nfeat;
}

int mdt_feature_angle_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type,
                          mdt_cb_feature_bond getbin, void *data,
                          mdt_cb_free freefunc)
{
  struct mdt_feature_bond *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_BOND, data, freefunc)->u.bond);
  feat->type = MDT_BOND_TYPE_ANGLE;
  feat->getbin = getbin;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_ANGLE, FALSE, 0);
  return nfeat;
}

int mdt_feature_dihedral_add(struct mdt_library *mlib, const char *name,
                             mod_mdt_calc precalc_type,
                             mdt_cb_feature_bond getbin, void *data,
                             mdt_cb_free freefunc)
{
  struct mdt_feature_bond *feat;
  int nfeat;

  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_BOND, data, freefunc)->u.bond);
  feat->type = MDT_BOND_TYPE_DIHEDRAL;
  feat->getbin = getbin;

  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_DIHEDRAL, FALSE, 0);
  return nfeat;
}

static gboolean check_feature(int ifeat, struct mdt_library *mlib, int nfeat,
                              GError **err)
{
  if (ifeat < 1 || ifeat > mlib->base.nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "Feature %d (%d) out of range 1 to %d", nfeat, ifeat,
                mlib->base.nfeat);
    return FALSE;
  } else {
    return TRUE;
  }
}

int mdt_feature_group_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, int ifeat1, int ifeat2,
                          mdt_cb_feature_group getbin, void *data,
                          mdt_cb_free freefunc, GError **err)
{
  struct mdt_feature_group *feat;
  int nfeat;

  if (!check_feature(ifeat1, mlib, 1, err)
      || !check_feature(ifeat2, mlib, 2, err)) {
    return -1;
  }
  feat = &(add_feature(mlib, &nfeat, MDT_FEATURE_GROUP, data,
                       freefunc)->u.group);
  feat->ifeat1 = ifeat1;
  feat->ifeat2 = ifeat2;
  feat->getbin = getbin;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_PROTEIN, FALSE, 0);
  return nfeat;
}
