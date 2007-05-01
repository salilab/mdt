/** \file mdt_tuples.h     Functions to build lists of atom tuples.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_TUPLE_H
#define __MDT_TUPLE_H

#include <glib.h>
#include "mdt_types.h"
#include "mod_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Get all tuples for a structure. */
G_GNUC_INTERNAL
struct mdt_tuple_list *tupclass(const struct structure *struc,
                                const struct sequence *seq,
                                const struct mdt_atom_class_list *atclass,
                                const struct libraries *libs);

G_END_DECLS

#endif  /* __MDT_TUPLE_H */
