/** \file mdt_tuples.h     Functions to build lists of atom tuples.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_TUPLE_H
#define __MDT_TUPLE_H

#include <glib.h>
#include "mdt_config.h"
#include "mdt_types.h"
#include "mod_types.h"

G_BEGIN_DECLS

/** Get all tuples for a structure. */
MDTDLLLOCAL
struct mdt_tuple_list *tupclass(const struct mod_structure *struc,
                                const struct mod_sequence *seq,
                                const struct mdt_atom_class_list *atclass,
                                const struct mod_libraries *libs);

G_END_DECLS

#endif  /* __MDT_TUPLE_H */
