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

/** A source of data for an MDT (generally an alignment) */
struct mdt_source;

/** Prepare a source alignment to add data to an MDT. Returns a source pointer
    (to be later freed with mdt_alignment_close()), or NULL on error. */
struct mdt_source *mdt_alignment_open(struct mod_mdt *mdt,
                                      const struct mdt_library *mlib,
                                      struct mod_alignment *aln, float distngh,
                                      gboolean sdchngh, int surftyp,
                                      int iacc1typ, struct mod_io_data *io,
                                      struct mod_libraries *libs, GError **err);

/** Close a source alignment previously opened with mdt_alignment_open(). */
void mdt_alignment_close(struct mdt_source *source);

/** Return the bin index (starting at 1) of a single MDT feature, at the
    given position in the source alignment. On failure, 0 is returned. */
int mdt_alignment_index(struct mdt_source *source, int ifeat, int is1, int ip1,
                        int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                        int ia1p, int ip2, int ibnd1, int ibnd1p, int is3,
                        int ir3, int ir3p, const struct mdt_library *mlib,
                        const struct mod_libraries *libs,
                        struct mod_energy_data *edat, GError **err);

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
