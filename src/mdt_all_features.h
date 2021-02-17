/** \file mdt_all_features.h    Functions to set up all types of feature.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
 */

#ifndef __MDT_ALL_FEATURES_H
#define __MDT_ALL_FEATURES_H

#include <glib.h>
#include <mod_types.h>
#include "mdt_config.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** Add a protein X-ray resolution feature. */
MDTDLLEXPORT
int mdt_feature_xray_resolution(struct mdt_library *mlib, int protein,
                                float nmr, GError **err);

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

/** Add an atom Z-coordinate feature. */
MDTDLLEXPORT
int mdt_feature_z_coordinate(struct mdt_library *mlib, gboolean pos2);

/** Add a fractional atom accessibility feature. */
MDTDLLEXPORT
int mdt_feature_fractional_atom_accessibility(struct mdt_library *mlib,
                                              gboolean pos2);

/** Add a Modeller atom type feature. */
MDTDLLEXPORT
int mdt_feature_atom_type(struct mdt_library *mlib, gboolean pos2);

/** Add a table atom feature. */
MDTDLLEXPORT
int mdt_feature_atom_table(struct mdt_library *mlib, gboolean pos2,
                           const char *table_name,
                           mdt_cb_get_property get_property,
                           gpointer data, GDestroyNotify freefunc);

/** Add an atom-atom distance feature. */
MDTDLLEXPORT
int mdt_feature_atom_distance(struct mdt_library *mlib);

/** Add an atom-atom bond separation feature. */
MDTDLLEXPORT
int mdt_feature_atom_bond_separation(struct mdt_library *mlib,
                                     gboolean disulfide);

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
                                         int protein, gboolean absolute,
                                         GError **err);

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

/** Add a residue-residue distance difference feature. */
MDTDLLEXPORT
int mdt_feature_residue_distance_difference(struct mdt_library *mlib,
                                            int protein1, int protein2,
                                            GError **err);

/** Add a mainchain conformation (Ramachandran) feature. */
MDTDLLEXPORT
int mdt_feature_mainchain_conformation(struct mdt_library *mlib, int protein,
                                       int delta, int align_delta,
                                       gboolean pos2,
                                       const struct mod_libraries *libs,
                                       GError **err);

/** Add a residue group feature. */
MDTDLLEXPORT
int mdt_feature_residue_group(struct mdt_library *mlib, int protein, int delta,
                              int align_delta, gboolean pos2,
                              int residue_grouping,
                              const struct mod_libraries *libs, GError **err);

/** Add a neighborhood difference feature. */
MDTDLLEXPORT
int mdt_feature_neighborhood_difference(struct mdt_library *mlib, int protein1,
                                        int protein2, GError **err);

/** Add an average neighborhood difference feature. */
MDTDLLEXPORT
int mdt_feature_average_neighborhood_difference(struct mdt_library *mlib,
                                                int protein1, int protein2,
                                                GError **err);

/** Add a sequence identity feature. */
MDTDLLEXPORT
int mdt_feature_sequence_identity(struct mdt_library *mlib, int protein1,
                                  int protein2, GError **err);

/** Add an average sidechain Biso feature. */
MDTDLLEXPORT
int mdt_feature_sidechain_biso(struct mdt_library *mlib, int protein,
                               int delta, int align_delta, gboolean pos2,
                               GError **err);

/** Add a protein alpha content feature. */
MDTDLLEXPORT
int mdt_feature_alpha_content(struct mdt_library *mlib, int protein,
                              GError **err);

/** Add a distance from a gap feature. */
MDTDLLEXPORT
int mdt_feature_gap_distance(struct mdt_library *mlib, int protein1,
                             int protein2, GError **err);

/** Add an average distance from a gap feature. */
MDTDLLEXPORT
int mdt_feature_average_gap_distance(struct mdt_library *mlib, int protein1,
                                     int protein2, GError **err);

/** Add a cluster feature. */
MDTDLLEXPORT
int mdt_feature_cluster(struct mdt_library *mlib, int ifeat1, int ifeat2,
                        int nbins, GError **err);

/** Add a bin pair to a cluster feature. Returns FALSE on error. */
MDTDLLEXPORT
gboolean mdt_cluster_add(struct mdt_library *mlib, int ifeat, int bin1,
                         int bin2, int bin, GError **err);

G_END_DECLS

#endif  /* __MDT_ALL_FEATURES_H */
