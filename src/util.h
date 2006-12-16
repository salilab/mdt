/** \file util.h           MDT utility functions.
 *
 *             Part of MODELLER, Copyright(c) 1989-2006 Andrej Sali
 */

#ifndef __MDT_UTIL_H
#define __MDT_UTIL_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
int *mdt_start_indices(const struct mdt_type *mdt);

/** Calculate the weights in the smoothing procedure for combining the
    a priori pdf with the experimental pdf. */
void weights(float weight, int nbins, float norm, float *w1, float *w2);

/** Return the index in the MDT pdf of the point whose feature bin
    indices are indf. */
int indmdt(const int *indf, const struct mdt_type *mdt);

/** Update the indices for the next point in the MDT. Return false if no
    more points are available. */
int roll_ind(int *indf, const struct mdt_type *mdt, int nfeat);

#ifdef __cplusplus
}
#endif
#endif  /* __MDT_UTIL_H */
