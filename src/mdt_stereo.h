/** \file mdt_stereo.h     Functions to determine stereochemistry for
 *                         template structures.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_STEREO_H
#define __MDT_STEREO_H

#include <glib.h>
#include "mdt_types.h"
#include "mod_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Get all of one type of bond (bond/angle/dihedral) for a structure. */
G_GNUC_INTERNAL
struct mdt_bond_list *get_stereo(const struct structure *struc,
                                 const struct sequence *seq,
                                 const struct mdt_atom_class_list *atclass,
                                 int bondtype, const struct libraries *libs);

G_END_DECLS

#endif  /* __MDT_STEREO_H */
