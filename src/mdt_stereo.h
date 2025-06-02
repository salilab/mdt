/** \file mdt_stereo.h     Functions to determine stereochemistry for
 *                         template structures.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#ifndef __MDT_STEREO_H
#define __MDT_STEREO_H

#include <glib.h>
#include "mdt_config.h"
#include "mdt_types.h"
#include "mod_types.h"

G_BEGIN_DECLS

/** Return True iff the atom index is OK, and the coordinates are
    defined. */
MDTDLLLOCAL
gboolean atmdefd(int ia1, const struct mod_coordinates *cd);

/** Get all of one type of bond (bond/angle/dihedral) for a structure. */
MDTDLLLOCAL
struct mdt_bond_list *get_stereo(const struct mod_structure *struc,
                                 const struct mod_sequence *seq,
                                 const struct mdt_atom_class_list *atclass,
                                 int bondtype,
                                 const struct mod_libraries *libs);

G_END_DECLS

#endif  /* __MDT_STEREO_H */
