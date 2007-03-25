/** \file mdt_index.h      Functions to calculate MDT indices.
 *
 *             Part of MODELLER, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_INDEX_H
#define __MDT_INDEX_H

#include <glib.h>
#include "mod_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Calculate a single MDT feature index */
G_GNUC_INTERNAL
int my_mdt_index(int ifi, const struct alignment *aln, int is1, int ip1,
                 int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                 int ia1p, const struct mdt_library *mlib, int ip2,
                 int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                 const struct libraries *libs,
                 const struct energy_data *edat, int *ierr);

G_END_DECLS

#endif  /* __MDT_INDEX_H */
