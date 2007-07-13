/** \file mdt_alignment.h      Functions to add alignment data to MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_ALIGNMENT_H
#define __MDT_ALIGNMENT_H

#include <glib.h>
#include "mdt_types.h"
#include "mod_types.h"

G_BEGIN_DECLS

/** Add data from an alignment to an MDT. Return TRUE on success. */
gboolean mdt_add_alignment(struct mod_mdt *mdt, const struct mdt_library *mlib,
                           struct mod_alignment *aln, float distngh,
                           gboolean sdchngh, int surftyp, int iacc1typ,
                           const int residue_span_range[4], int pairs,
                           int triples, struct mod_io_data *io,
                           struct mod_energy_data *edat,
                           struct mod_libraries *libs, GError **err);

G_END_DECLS

#endif  /* __MDT_ALIGNMENT_H */
