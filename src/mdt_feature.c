/** \file mdt_feature.c    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "mdt_feature.h"
#include "mdt_index.h"

/** Is the given feature type periodic? */
gboolean mdt_feature_is_periodic(int ifeat)
{
  switch (ifeat) {
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
  case 28:
  case 29:
  case 41:
  case 42:
  case 53:
  case 54:
  case 55:
  case 56:
  case 57:
  case 106:
  case 107:
  case 108:
  case 114:
    return TRUE;
  default:
    return FALSE;
  }
}

void mdt_feature_add_needed_file(struct mdt_library *mlib, int ifeat,
                                 mod_mdt_file filetype)
{
  struct mod_mdt_libfeature *feat = &mlib->base.features[ifeat - 1];
  feat->ndatfil++;
  feat->idatfil = g_realloc(feat->idatfil, sizeof(int) * feat->ndatfil);
  feat->idatfil[feat->ndatfil - 1] = filetype;
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

/** \return TRUE iff 'protein' is in range. */
static gboolean check_protein(int protein, const char *feattype, GError **err)
{
  if (protein < 0 || protein > 1) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s features can act only on protein 0 or protein 1; "
                "%d was given", feattype, protein);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Helper function to add a new feature structure, and return it. */
static struct mdt_feature *add_feature(struct mdt_library *mlib, int *nfeat)
{
  mlib->feature_added = TRUE;
  *nfeat = mlib->base.nfeat + 1;
  mlib->features = g_array_set_size(mlib->features, *nfeat);
  return &g_array_index(mlib->features, struct mdt_feature, *nfeat - 1);
}

int mdt_feature_protein_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein,
                            mdt_cb_feature_protein getbin, void *data,
                            GError **err)
{
  GString *fullname;
  struct mdt_feature *feat;
  int nfeat;

  if (!check_protein(protein, "Protein", err)) {
    return -1;
  }

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_PROTEIN;
  feat->u.protein.protein = protein;
  feat->u.protein.getbin = getbin;
  feat->data = data;
  fullname = g_string_new(name);
  g_string_append_printf(fullname, " of protein %d", protein);
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              protein == 0 ? MOD_MDTP_A : MOD_MDTP_B,
                              MOD_MDTS_PROTEIN, FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}

int mdt_feature_residue_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein, int delta,
                            int align_delta, gboolean pos2,
                            int bin_seq_outrange, mdt_cb_feature_residue getbin,
                            void *data, GError **err)
{
  GString *fullname;
  struct mdt_feature *feat;
  int nfeat;

  if (!check_protein(protein, "Residue", err)) {
    return -1;
  }

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_RESIDUE;
  feat->u.residue.protein = protein;
  feat->u.residue.delta = delta;
  feat->u.residue.align_delta = align_delta;
  feat->u.residue.pos2 = pos2;
  feat->u.residue.bin_seq_outrange = bin_seq_outrange;
  feat->u.residue.getbin = getbin;
  feat->data = data;
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
                              protein == 0 ? MOD_MDTP_A : MOD_MDTP_B,
                              pos2 ? MOD_MDTS_RESIDUE_PAIR : MOD_MDTS_RESIDUE,
                              FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}

int mdt_feature_residue_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein,
                                 gboolean asymmetric,
                                 mdt_cb_feature_residue_pair getbin, void *data,
                                 GError **err)
{
  struct mdt_feature *feat;
  int nfeat;

  if (!check_protein(protein, "Residue pair", err)) {
    return -1;
  }

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_RESIDUE_PAIR;
  feat->u.residue_pair.protein = protein;
  feat->u.residue_pair.getbin = getbin;
  feat->data = data;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              protein == 0 ? MOD_MDTP_A : MOD_MDTP_B,
                              MOD_MDTS_RESIDUE_PAIR, asymmetric, 0);
  return nfeat;
}

int mdt_feature_atom_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type, gboolean pos2,
                         mdt_cb_feature_atom getbin, void *data)
{
  GString *fullname;
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_ATOM;
  feat->u.atom.pos2 = pos2;
  feat->u.atom.getbin = getbin;
  feat->data = data;
  fullname = g_string_new(name);
  if (pos2) {
    g_string_append(fullname, " at pos2");
  }
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              MOD_MDTP_A,
                              pos2 ? MOD_MDTS_ATOM_PAIR : MOD_MDTS_ATOM,
                              FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}

int mdt_feature_atom_pair_add(struct mdt_library *mlib, const char *name,
                              mod_mdt_calc precalc_type, gboolean asymmetric,
                              mdt_cb_feature_atom_pair getbin, void *data)
{
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_ATOM_PAIR;
  feat->u.atom_pair.getbin = getbin;
  feat->data = data;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_ATOM_PAIR, asymmetric, 0);
  return nfeat;
}

int mdt_feature_tuple_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, gboolean pos2,
                          mdt_cb_feature_tuple getbin, void *data)
{
  GString *fullname;
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_TUPLE;
  feat->u.tuple.pos2 = pos2;
  feat->u.tuple.getbin = getbin;
  feat->data = data;
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
                               mdt_cb_feature_tuple_pair getbin, void *data)
{
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_TUPLE_PAIR;
  feat->u.tuple_pair.getbin = getbin;
  feat->data = data;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR, TRUE, 0);
  return nfeat;
}

int mdt_feature_bond_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type,
                         mdt_cb_feature_bond getbin, void *data)
{
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_BOND;
  feat->u.bond.type = MDT_BOND_TYPE_BOND;
  feat->u.bond.getbin = getbin;
  feat->data = data;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_BOND, FALSE, 0);
  return nfeat;
}

int mdt_feature_angle_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type,
                          mdt_cb_feature_bond getbin, void *data)
{
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_BOND;
  feat->u.bond.type = MDT_BOND_TYPE_ANGLE;
  feat->u.bond.getbin = getbin;
  feat->data = data;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_ANGLE, FALSE, 0);
  return nfeat;
}

int mdt_feature_dihedral_add(struct mdt_library *mlib, const char *name,
                             mod_mdt_calc precalc_type,
                             mdt_cb_feature_bond getbin, void *data)
{
  struct mdt_feature *feat;
  int nfeat;

  feat = add_feature(mlib, &nfeat);
  feat->type = MDT_FEATURE_BOND;
  feat->u.bond.type = MDT_BOND_TYPE_DIHEDRAL;
  feat->u.bond.getbin = getbin;
  feat->data = data;
  mod_mdt_libfeature_register(&mlib->base, nfeat, name, precalc_type,
                              MOD_MDTP_A, MOD_MDTS_DIHEDRAL, FALSE, 0);
  return nfeat;
}
