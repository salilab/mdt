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

/** Add a residue feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_residue_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein, int delta,
                            gboolean pos2, mdt_cb_feature_residue getbin,
                            void *data, GError **err);

/** Add an atom feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_atom_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type, gboolean pos2,
                         mdt_cb_feature_atom getbin, void *data);

/** Add a tuple feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_tuple_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, gboolean pos2,
                          mdt_cb_feature_tuple getbin, void *data);

/** Add a data file type needed for a given feature */
MDTDLLEXPORT
void mdt_feature_add_needed_file(struct mdt_library *mlib, int ifeat,
                                 mod_mdt_file filetype);

/** Set the number of bins for a given feature.
    \note One extra bin is always created - the 'undefined' bin. */
MDTDLLEXPORT
void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins);

/** Set the range and symbol for a given bin and feature. */
MDTDLLEXPORT
void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol);

/** Add a protein X-ray resolution feature. */
MDTDLLEXPORT
int mdt_feature_xray_resolution(struct mdt_library *mlib, int protein,
                                GError **err);

/** Add a protein radius of gyration feature. */
MDTDLLEXPORT
int mdt_feature_radius_of_gyration(struct mdt_library *mlib, int protein,
                                   GError **err);

/** Add a residue accessibility feature. */
MDTDLLEXPORT
int mdt_feature_residue_accessibility(struct mdt_library *mlib, int protein,
                                      int delta, gboolean pos2, GError **err);

/** Add an atom accessibility feature. */
MDTDLLEXPORT
int mdt_feature_atom_accessibility(struct mdt_library *mlib, gboolean pos2);

/** Add a fractional atom accessibility feature. */
MDTDLLEXPORT
int mdt_feature_fractional_atom_accessibility(struct mdt_library *mlib,
                                              gboolean pos2);

/** Add a Modeller atom type feature. */
MDTDLLEXPORT
int mdt_feature_atom_type(struct mdt_library *mlib, gboolean pos2);

/** Add a hydrogen bond donor feature. */
MDTDLLEXPORT
int mdt_feature_hydrogen_bond_donor(struct mdt_library *mlib, gboolean pos2);

/** Add a hydrogen bond acceptor feature. */
MDTDLLEXPORT
int mdt_feature_hydrogen_bond_acceptor(struct mdt_library *mlib, gboolean pos2);

/** Add a total hydrogen bond charge around an atom feature. */
MDTDLLEXPORT
int mdt_feature_hydrogen_bond_charge(struct mdt_library *mlib, gboolean pos2);

/** Add a protein hydrogen bond satisfaction feature. */
MDTDLLEXPORT
int mdt_feature_hydrogen_bond_satisfaction(struct mdt_library *mlib,
                                           int protein, GError **err);

/** Add a tuple type feature. */
MDTDLLEXPORT
int mdt_feature_tuple_type(struct mdt_library *mlib, gboolean pos2);

G_END_DECLS

#endif  /* __MDT_FEATURE_H */
