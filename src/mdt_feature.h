/** \file mdt_feature.h    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_FEATURE_H
#define __MDT_FEATURE_H

#include <glib.h>
#include <mod_types.h>
#include "mdt_config.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** Add a protein feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_protein_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein,
                            mdt_cb_feature_protein getbin, void *data,
                            GError **err);

/** Add a protein pair feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_protein_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein1,
                                 int protein2,
                                 mdt_cb_feature_protein_pair getbin, void *data,
                                 GError **err);

/** Add a residue feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_residue_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein, int delta,
                            int align_delta, gboolean pos2,
                            int bin_seq_outrange, mdt_cb_feature_residue getbin,
                            void *data, GError **err);

/** Add a residue pair feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_residue_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein,
                                 gboolean asymmetric,
                                 mdt_cb_feature_residue_pair getbin, void *data,
                                 GError **err);

/** Add an aligned residue feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_aligned_residue_add(struct mdt_library *mlib, const char *name,
                                    mod_mdt_calc precalc_type, int protein1,
                                    int protein2,
                                    mdt_cb_feature_aligned_residue getbin,
                                    void *data, GError **err);

/** Add an aligned residue pair feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_aligned_residue_pair_add(
    struct mdt_library *mlib, const char *name, mod_mdt_calc precalc_type,
    int protein1, int protein2, gboolean asymmetric,
    mdt_cb_feature_aligned_residue_pair getbin, void *data, GError **err);

/** Add an atom feature.
    \note The system is automatically instructed to read in PDB files.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_atom_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type, gboolean pos2,
                         mdt_cb_feature_atom getbin, void *data);

/** Add an atom pair feature.
    \note The system is automatically instructed to read in PDB files.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_atom_pair_add(struct mdt_library *mlib, const char *name,
                              mod_mdt_calc precalc_type, gboolean asymmetric,
                              mdt_cb_feature_atom_pair getbin, void *data);

/** Add a tuple feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_tuple_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, gboolean pos2,
                          mdt_cb_feature_tuple getbin, void *data);

/** Add a tuple pair feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_tuple_pair_add(struct mdt_library *mlib, const char *name,
                               mod_mdt_calc precalc_type,
                               mdt_cb_feature_tuple_pair getbin, void *data);

/** Add a bond feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_bond_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type,
                         mdt_cb_feature_bond getbin, void *data);

/** Add an angle feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_angle_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type,
                          mdt_cb_feature_bond getbin, void *data);

/** Add a dihedral feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_dihedral_add(struct mdt_library *mlib, const char *name,
                             mod_mdt_calc precalc_type,
                             mdt_cb_feature_bond getbin, void *data);

/** Add a data file type needed for a given feature */
MDTDLLEXPORT
void mdt_feature_add_needed_file(struct mdt_library *mlib, int ifeat,
                                 mod_mdt_file filetype);

/** Set whether a feature's range is periodic (FALSE by default) */
MDTDLLEXPORT
void mdt_feature_periodic_set(struct mdt_library *mlib, int ifeat,
                              gboolean periodic);

/** Get whether a feature's range is periodic */
MDTDLLEXPORT
gboolean mdt_feature_periodic_get(const struct mdt_library *mlib, int ifeat);

/** Set the number of bins for a given feature.
    \note One extra bin is always created - the 'undefined' bin. */
MDTDLLEXPORT
void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins);

/** Set the range and symbol for a given bin and feature. */
MDTDLLEXPORT
void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol);

G_END_DECLS

#endif  /* __MDT_FEATURE_H */
