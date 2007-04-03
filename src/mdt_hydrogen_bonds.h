/** \file mdt_hydrogen_bonds.h    Functions for handling hydrogen bonds
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_HBOND_H
#define __MDT_HBOND_H

#include <glib.h>
#include "mdt_types.h"
#include "mod_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
G_GNUC_INTERNAL
int numb_hda(int ia, const int hb_iattyp[], const struct coordinates *cd,
             const struct mdt_atom_class_list *atclass, float hbond_cutoff,
             gboolean acceptor, int nbins);

G_END_DECLS

#endif  /* __MDT_HBOND_H */
