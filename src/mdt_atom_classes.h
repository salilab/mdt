/** \file mdt_atom_classes.h  Functions to handle atom classes.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_ATOMCLASS_H
#define __MDT_ATOMCLASS_H

#include <glib.h>
#include "mdt_types.h"

G_BEGIN_DECLS

/** Make a new atom class list. */
struct mdt_atom_class_list *mdt_atom_class_list_new(int natom);

/** Free an existing atom class list. */
void mdt_atom_class_list_free(struct mdt_atom_class_list *atclass);

/** Read atom class information from a file; return TRUE on success. */
gboolean mdt_atom_classes_read(const gchar *filename,
                               struct mdt_library *mlib, int natom,
                               GError **err);

/** Read hydrogen bond class information from a file; return TRUE on success. */
gboolean mdt_hbond_read(const gchar *filename, struct mdt_library *mlib, 
                        GError **err);

G_END_DECLS

#endif  /* __MDT_ATOMCLASS_H */
