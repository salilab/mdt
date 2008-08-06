/** \file mdt_feature.h    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_FEATURE_H
#define __MDT_FEATURE_H

#include <glib.h>
#include <mod_types.h>
#include "mdt_config.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** Add a protein feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_protein_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein,
                            mdt_cb_feature_protein getbin, void *data,
                            GError **err);

/** Set the number of bins for a given feature.
    \note One extra bin is always created - the 'undefined' bin. */
MDTDLLEXPORT
void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins);

/** Set the range and symbol for a given bin and feature. */
MDTDLLEXPORT
void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol);

/** Add a protein X-ray resolution feature. */
MDTDLLEXPORT
int mdt_feature_xray_resolution(struct mdt_library *mlib, int protein,
                                GError **err);

/** Add a protein radius of gyration feature. */
MDTDLLEXPORT
int mdt_feature_radius_of_gyration(struct mdt_library *mlib, int protein,
                                   GError **err);

G_END_DECLS

#endif  /* __MDT_FEATURE_H */
