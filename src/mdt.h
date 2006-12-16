/** \file mdt.h            Functions to handle MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#ifndef __MDT_MDT_H
#define __MDT_MDT_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Smooth a histogram or the 2D plot with a uniform prior */
void mdt_smooth(const struct mdt_type *mdtin, struct mdt_type *mdtout,
                int dimensions, float weight, int *ierr);

/** Transform an MDT with an exp function */
void mdt_exp_transform(struct mdt_type *mdt, float offset, float expoffset,
                       float multiplier, float power);

/** Transform an MDT with a linear function */
void mdt_linear_transform(struct mdt_type *mdt, float offset, float multiplier);

/** Transform an MDT with an inverse function */
void mdt_inverse_transform(struct mdt_type *mdt, float offset, float multiplier,
                           float undefined);

#ifdef __cplusplus
}
#endif
#endif  /* __MDT_MDT_H */
