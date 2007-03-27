/** \file mdt_types.h      Functions to handle MDT types.
 *
 *             Part of MODELLER, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_TYPES_H
#define __MDT_TYPES_H

#include <glib.h>
#include "mod_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Library of feature data used by MDTs */
struct mdt_library {
  /** Base Modeller type */
  struct mod_mdt_library base;
  /** Deltas for some feature types ('residue + delta i') */
  int deltai, deltaj;
  /** TRUE if deltas refer to align. positions, or FALSE if residue positions */
  gboolean deltai_ali, deltaj_ali;
};

/** Make a new mdt_library structure */
struct mdt_library *mdt_library_new(void);

/** Free an mdt_library structure */
void mdt_library_free(struct mdt_library *mlib);

G_END_DECLS

#endif  /* __MDT_TYPES_H */
