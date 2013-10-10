/** \file mdt_property.h  Functions to precalculate protein properties used
 *                        to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#ifndef __MDT_PROPERTY_H
#define __MDT_PROPERTY_H

#include <glib.h>
#include "mdt_config.h"
#include "mod_types.h"
#include "mdt_index.h"

G_BEGIN_DECLS

/** A single disulfide bridge */
struct mdt_disulfide {
  /** Indices of the two atoms */
  int atom1, atom2;
  /** Indices of the two residues */
  int res1, res2;
};

/** A list of disulfide bridges */
struct mdt_disulfide_list {
  /** Number of bridges */
  int nss;
  /** Atom and residue information for each bridge */
  struct mdt_disulfide *ss;
};

/** Precalculated per-sequence properties for calculating MDT indices.
    Note that this is an 'opaque' struct; it should only be accessed
    via the functions in this header. */
struct mdt_properties;

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

/** Get/calculate all excluded atom pairs for a structure.
    NULL is returned if no pairs are excluded. */
MDTDLLLOCAL
gboolean property_exclusions(const struct mod_alignment *aln,
                             int is, struct mdt_properties *prop,
                             const struct mdt_library *mlib,
                             gboolean exclude_bonds,
                             gboolean exclude_angles,
                             gboolean exclude_dihedrals,
                             const struct mod_libraries *libs,
                             GHashTable **exclusions, GError **err);

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
                             GError **err);

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
                          const struct mdt_library *mlib,
                          const struct mod_libraries *libs, GError **err);

/** Get/calculate the array of hydrogen bond atom type bin indices */
MDTDLLLOCAL
const int *property_hb_iatta(const struct mod_alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib,
                             const struct mod_libraries *libs, GError **err);

/** Get/calculate the hydrogen bond satisfaction index */
MDTDLLLOCAL
gboolean property_hbpot(const struct mod_alignment *aln, int is,
                        struct mdt_properties *prop,
                        const struct mdt_library *mlib,
                        const struct mod_libraries *libs, float *hbpot,
                        GError **err);

/** Get/calculate the array of fractional atom accessibilities.
    \return NULL on failure. */
MDTDLLLOCAL
const float *property_fatmacc(const struct mod_alignment *aln, int is,
                              struct mdt_properties *prop,
                              const struct mod_libraries *libs, GError **err);

/** Get/calculate the radius of gyration */
MDTDLLLOCAL
float property_radius_gyration(const struct mod_alignment *aln, int is,
                               struct mdt_properties *prop);

/** Get/calculate the array of average sidechain Biso. */
MDTDLLLOCAL
const float *property_sidechain_biso(const struct mod_alignment *aln, int is,
                                     struct mdt_properties *prop);

/** Get/calculate the array of per-residue distance atom indices */
MDTDLLLOCAL
void property_distance_atom_indices(const struct mod_alignment *aln, int is,
                                    struct mdt_properties *prop,
                                    const struct mdt_library *mlib,
                                    const int **dstind1, const int **dstind2);

/** Get/calculate the array of atom types for residue bond separation */
MDTDLLLOCAL
const int *property_resbond_attyp(const struct mod_alignment *aln, int is,
                                  struct mdt_properties *prop,
                                  const struct mdt_library *mlib,
                                  const struct mod_libraries *libs);

/** Get/calculate the list of disulfide bridges */
MDTDLLLOCAL
const struct mdt_disulfide_list *property_disulfides(
                               const struct mod_alignment *aln, int is,
                               struct mdt_properties *prop,
                               const struct mdt_library *mlib,
                               const struct mod_libraries *libs);

/** Copy any write callbacks to an MDT table */
MDTDLLLOCAL
void mdt_property_get_write_callbacks(const struct mdt_properties *prop,
                                      const struct mod_alignment *aln,
                                      struct mdt *mdt);

/* Pack two atom or bin indices into a pointer.
   Note  - no implementation for platforms where pointer is not 32 or 64 bits.
         - indices are first made unsigned, then padded out to the size
           of the pointer (two casts).
         - on 32-bit platforms, indices bigger than 2^16 cannot be stored;
           a runtime error is emitted on such a system.
 */
#if MDT_SIZEOF_POINTER == 8
#define MAKE_HASH_KEY_ASYMMETRIC(a, b) \
        ((gpointer)((guint64)(guint32)(a) << 32 | (guint64)(guint32)(b)))
#define GET_HASH_KEY_LOW(h) \
        ((guint32)(guint64)h)
#define GET_HASH_KEY_HIGH(h) \
        ((guint32)((guint64)h >> 32))
#define MAKE_HASH_KEY(a, b) \
        (a < b ? (gpointer)((guint64)(guint32)a << 32 | (guint64)(guint32)b) \
               : (gpointer)((guint64)(guint32)b << 32 | (guint64)(guint32)a))
#elif MDT_SIZEOF_POINTER == 4
#define MAKE_HASH_KEY_ASYMMETRIC(a, b) \
        ((gpointer)((guint32)(guint16)(a) << 16 | (guint32)(guint16)(b)))
#define GET_HASH_KEY_LOW(h) \
        ((guint16)(guint32)h)
#define GET_HASH_KEY_HIGH(h) \
        ((guint16)((guint32)h >> 16))
#define MAKE_HASH_KEY(a, b) \
        (a < b ? (gpointer)((guint32)(guint16)a << 16 | (guint32)(guint16)b) \
               : (gpointer)((guint32)(guint16)b << 16 | (guint32)(guint16)a))
#endif

G_END_DECLS

#endif  /* __MDT_PROPERTY_H */
