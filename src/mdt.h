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

/** Transform an MDT with a log function */
void mdt_log_transform(struct mdt_type *mdt, float offset, float multiplier,
                       float undefined);

/** Offset an MDT by the minimum value. */
void mdt_offset_min(struct mdt_type *mdt, int dimensions, int *ierr);

/** Close an MDT so that it is useful for creating periodic splines. */
void mdt_close(struct mdt_type *mdt, int dimensions, int *ierr);

/** Get the entropy of the dependent variable. */
float mdt_entropy_hx(const struct mdt_type *mdt, int *ierr);

/** Write out an MDT. */
void mdt_write(const struct mdt_type *mdt, const struct mdt_library *mlib,
               const char *filename, mbool write_preamble, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* __MDT_MDT_H */
