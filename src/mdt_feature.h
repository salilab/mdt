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

/** Add a protein sequence length feature. */
MDTDLLEXPORT
int mdt_feature_sequence_length(struct mdt_library *mlib, int protein,
                                GError **err);

/** Add a residue accessibility feature. */
MDTDLLEXPORT
int mdt_feature_residue_accessibility(struct mdt_library *mlib, int protein,
                                      int delta, int align_delta, gboolean pos2,
                                      GError **err);

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

/** Add an atom-atom distance feature. */
MDTDLLEXPORT
int mdt_feature_atom_distance(struct mdt_library *mlib);

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

/** Add a tuple-tuple non-bonded distance feature. */
MDTDLLEXPORT
int mdt_feature_tuple_distance(struct mdt_library *mlib);

/** Add a tuple-tuple non-bonded angle1 feature. */
MDTDLLEXPORT
int mdt_feature_tuple_angle1(struct mdt_library *mlib, GError **err);

/** Add a tuple-tuple non-bonded angle2 feature. */
MDTDLLEXPORT
int mdt_feature_tuple_angle2(struct mdt_library *mlib, GError **err);

/** Add a tuple-tuple non-bonded dihedral1 feature. */
MDTDLLEXPORT
int mdt_feature_tuple_dihedral1(struct mdt_library *mlib, GError **err);

/** Add a tuple-tuple non-bonded dihedral2 feature. */
MDTDLLEXPORT
int mdt_feature_tuple_dihedral2(struct mdt_library *mlib, GError **err);

/** Add a tuple-tuple non-bonded dihedral3 feature. */
MDTDLLEXPORT
int mdt_feature_tuple_dihedral3(struct mdt_library *mlib, GError **err);

/** Add a bond type feature. */
MDTDLLEXPORT
int mdt_feature_bond_type(struct mdt_library *mlib);

/** Add an angle type feature. */
MDTDLLEXPORT
int mdt_feature_angle_type(struct mdt_library *mlib);

/** Add a dihedral type feature. */
MDTDLLEXPORT
int mdt_feature_dihedral_type(struct mdt_library *mlib);

/** Add a bond length feature. */
MDTDLLEXPORT
int mdt_feature_bond_length(struct mdt_library *mlib);

/** Add an angle feature. */
MDTDLLEXPORT
int mdt_feature_angle(struct mdt_library *mlib);

/** Add a dihedral feature. */
MDTDLLEXPORT
int mdt_feature_dihedral(struct mdt_library *mlib);

/** Add a chi1 dihedral feature. */
MDTDLLEXPORT
int mdt_feature_chi1_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err);

/** Add a chi2 dihedral feature. */
MDTDLLEXPORT
int mdt_feature_chi2_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err);

/** Add a chi3 dihedral feature. */
MDTDLLEXPORT
int mdt_feature_chi3_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err);

/** Add a chi4 dihedral feature. */
MDTDLLEXPORT
int mdt_feature_chi4_dihedral(struct mdt_library *mlib, int protein,
                              int delta, int align_delta, gboolean pos2,
                              GError **err);

/** Add a phi dihedral feature. */
MDTDLLEXPORT
int mdt_feature_phi_dihedral(struct mdt_library *mlib, int protein,
                             int delta, int align_delta, gboolean pos2,
                             GError **err);

/** Add a psi dihedral feature. */
MDTDLLEXPORT
int mdt_feature_psi_dihedral(struct mdt_library *mlib, int protein,
                             int delta, int align_delta, gboolean pos2,
                             GError **err);

/** Add a omega dihedral feature. */
MDTDLLEXPORT
int mdt_feature_omega_dihedral(struct mdt_library *mlib, int protein,
                               int delta, int align_delta, gboolean pos2,
                               GError **err);

/** Add an alpha dihedral feature. */
MDTDLLEXPORT
int mdt_feature_alpha_dihedral(struct mdt_library *mlib, int protein,
                               int delta, int align_delta, gboolean pos2,
                               GError **err);

/** Add a chi1 dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_chi1_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err);

/** Add a chi2 dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_chi2_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err);

/** Add a chi3 dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_chi3_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err);

/** Add a chi4 dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_chi4_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err);

/** Add a chi5 dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_chi5_class(struct mdt_library *mlib, int protein,
                           int delta, int align_delta, gboolean pos2,
                           const struct mod_libraries *libs, GError **err);

/** Add a phi dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_phi_class(struct mdt_library *mlib, int protein,
                          int delta, int align_delta, gboolean pos2,
                          const struct mod_libraries *libs, GError **err);

/** Add a psi dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_psi_class(struct mdt_library *mlib, int protein,
                          int delta, int align_delta, gboolean pos2,
                          const struct mod_libraries *libs, GError **err);

/** Add an omega dihedral class feature. */
MDTDLLEXPORT
int mdt_feature_omega_class(struct mdt_library *mlib, int protein,
                            int delta, int align_delta, gboolean pos2,
                            const struct mod_libraries *libs, GError **err);

/** Add a residue type feature. */
MDTDLLEXPORT
int mdt_feature_residue_type(struct mdt_library *mlib, int protein,
                             int delta, int align_delta, gboolean pos2,
                             const struct mod_libraries *libs, GError **err);

/** Add a residue-residue distance feature. */
MDTDLLEXPORT
int mdt_feature_residue_distance(struct mdt_library *mlib, int protein,
                                 GError **err);

/** Add an average accessibility of a residue pair feature. */
MDTDLLEXPORT
int mdt_feature_average_residue_accessibility(struct mdt_library *mlib,
                                              int protein, GError **err);

/** Add a residue index difference feature. */
MDTDLLEXPORT
int mdt_feature_residue_index_difference(struct mdt_library *mlib,
                                         int protein, GError **err);

/** Add a psi dihedral difference feature. */
MDTDLLEXPORT
int mdt_feature_psi_dihedral_difference(struct mdt_library *mlib,
                                        int protein1, int protein2,
                                        GError **err);

/** Add a phi dihedral difference feature. */
MDTDLLEXPORT
int mdt_feature_phi_dihedral_difference(struct mdt_library *mlib,
                                        int protein1, int protein2,
                                        GError **err);

/** Add an omega dihedral difference feature. */
MDTDLLEXPORT
int mdt_feature_omega_dihedral_difference(struct mdt_library *mlib,
                                          int protein1, int protein2,
                                          GError **err);

G_END_DECLS

#endif  /* __MDT_FEATURE_H */
