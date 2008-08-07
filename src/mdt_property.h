/** \file mdt_property.h  Functions to precalculate protein properties used
 *                        to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_PROPERTY_H
#define __MDT_PROPERTY_H

#include <glib.h>
#include "mdt_config.h"
#include "mod_types.h"
#include "mdt_index.h"

G_BEGIN_DECLS

/** Make a new mdt_properties structure */
MDTDLLLOCAL
struct mdt_properties *mdt_properties_new(const struct mod_alignment *aln);

/** Free an mdt_properties structure */
MDTDLLLOCAL
void mdt_properties_free(struct mdt_properties *prop,
                         const struct mod_alignment *aln);

/** Get/calculate the list of all bonds for a structure. */
MDTDLLLOCAL
const struct mdt_bond_list *property_bonds(const struct mod_alignment *aln,
                                           int is,
                                           struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int bondtype,
                                           const struct mod_libraries *libs);

/** Get a single bond from a structure */
MDTDLLLOCAL
const struct mdt_bond *property_one_bond(const struct mod_alignment *aln,
                                         int is, struct mdt_properties *prop,
                                         const struct mdt_library *mlib,
                                         int bondtype, int ibnd1,
                                         const struct mod_libraries *libs);

/** Get/calculate the list of all tuples for a structure. */
MDTDLLLOCAL
const struct mdt_tuple_list *property_tuples(const struct mod_alignment *aln,
                                             int is,
                                             struct mdt_properties *prop,
                                             const struct mdt_library *mlib,
                                             const struct mod_libraries *libs);

/** Require that the tuples have at least min_natom atoms each */
MDTDLLLOCAL
gboolean tuple_require_natom(const struct mdt_library *mlib, int min_natom,
                             int ifeat, GError **err);

/** Get a single atom tuple from a structure */
MDTDLLLOCAL
const struct mdt_tuple *property_one_tuple(const struct mod_alignment *aln,
                                           int is,
                                           struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int ibnd1, int ia1,
                                           const struct mod_libraries *libs);

/** Get/calculate the array of atom type bin indices */
MDTDLLLOCAL
const int *property_iatta(const struct mod_alignment *aln, int is,
                          struct mdt_properties *prop,
                          const struct mdt_library *mlib, int ifi,
                          const struct mod_libraries *libs, GError **err);

/** Get/calculate the array of hydrogen bond atom type bin indices */
MDTDLLLOCAL
const int *property_hb_iatta(const struct mod_alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct mod_libraries *libs, GError **err);

/** Get/calculate the hydrogen bond satisfaction index */
MDTDLLLOCAL
gboolean property_hbpot(const struct mod_alignment *aln, int is,
                        struct mdt_properties *prop,
                        const struct mdt_library *mlib, int ifi,
                        const struct mod_libraries *libs, float *hbpot,
                        GError **err);

/** Get/calculate the array of fractional atom accessibilities.
    \return TRUE on success. */
MDTDLLLOCAL
gboolean property_fatmacc(const struct mod_alignment *aln, int is,
                          struct mdt_properties *prop,
                          const struct mod_mdt_libfeature *feat,
                          const struct mod_libraries *libs, float **table,
                          GError **err);

/** Get/calculate the radius of gyration */
MDTDLLLOCAL
float property_radius_gyration(const struct mod_alignment *aln, int is,
                               struct mdt_properties *prop);

G_END_DECLS

#endif  /* __MDT_PROPERTY_H */
