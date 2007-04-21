/** \file mdt_property.h  Functions to precalculate protein properties used
 *                        to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_PROPERTY_H
#define __MDT_PROPERTY_H

#include <glib.h>
#include "mod_types.h"
#include "mdt_index.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Make a new mdt_properties structure */
G_GNUC_INTERNAL
struct mdt_properties *mdt_properties_new(const struct alignment *aln);

/** Free an mdt_properties structure */
G_GNUC_INTERNAL
void mdt_properties_free(struct mdt_properties *prop,
                         const struct alignment *aln);

/** Get/calculate the list of all bonds for a structure. */
G_GNUC_INTERNAL
const struct mdt_bond_list *property_bonds(const struct alignment *aln, int is,
                                           struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int bondtype,
                                           const struct libraries *libs);

/** Get a single bond from a structure */
G_GNUC_INTERNAL
const struct mdt_bond *property_one_bond(const struct alignment *aln,
                                         int is, struct mdt_properties *prop,
                                         const struct mdt_library *mlib,
                                         int bondtype, int ibnd1,
                                         const struct libraries *libs);

/** Get/calculate the list of all triplets for a structure. */
G_GNUC_INTERNAL
const struct mdt_triplet_list *property_triplets(const struct alignment *aln,
                                                 int is,
                                                 struct mdt_properties *prop,
                                                 const struct mdt_library *mlib,
                                                 const struct libraries *libs);

/** Get a single atom triplet from a structure */
G_GNUC_INTERNAL
const struct mdt_triplet *property_one_triplet(const struct alignment *aln,
                                               int is,
                                               struct mdt_properties *prop,
                                               const struct mdt_library *mlib,
                                               int ibnd1, int ia1,
                                               const struct libraries *libs);

/** Get/calculate the resolution bin index */
G_GNUC_INTERNAL
int property_iresol(const struct alignment *aln, int is,
                    struct mdt_properties *prop,
                    const struct mdt_library *mlib, int ifi,
                    const struct mdt_libfeature *feat);

/** Get/calculate the array of atom type bin indices */
G_GNUC_INTERNAL
const int *property_iatta(const struct alignment *aln, int is,
                          struct mdt_properties *prop,
                          const struct mdt_library *mlib, int ifi,
                          const struct libraries *libs, GError **err);

/** Get/calculate the array of hydrogen bond atom type bin indices */
G_GNUC_INTERNAL
const int *property_hb_iatta(const struct alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct libraries *libs, GError **err);

/** Get/calculate the hydrogen bond satisfaction index */
G_GNUC_INTERNAL
gboolean property_hbpot(const struct alignment *aln, int is,
                        struct mdt_properties *prop,
                        const struct mdt_library *mlib, int ifi,
                        const struct libraries *libs, float *hbpot,
                        GError **err);

/** Get/calculate the array of atom accessibility bin indices */
G_GNUC_INTERNAL
const int *property_iatmacc(const struct alignment *aln, int is,
                            struct mdt_properties *prop,
                            const struct mdt_library *mlib, int ifi,
                            const struct mdt_libfeature *feat);

/** Get/calculate the array of fractional atom accessibility bin indices */
G_GNUC_INTERNAL
const int *property_ifatmacc(const struct alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct mdt_libfeature *feat,
                             const struct libraries *libs, GError **err);

/** Get/calculate the radius of gyration bin index */
G_GNUC_INTERNAL
int property_radius_gyration(const struct alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct mdt_libfeature *feat);

G_END_DECLS

#endif  /* __MDT_PROPERTY_H */
