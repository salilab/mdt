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
int roll_ind(int indf[], const int istart[], const int iend[], int nfeat);

/** Roll n indices in ind so that all combinations of n different indices
    are generated, where the possible values for each index are 1+x to nmax-y.

    For example, for n=3:

    for (i = 0, i < nmax-2; i++) {
      for (j = i+1, j < nmax-1; j++) {
        for (k = j+1, k < nmax; k++) {

    *ind should be NULL on the first call to this routine, and it will be
    initialized (the user should free it when finished).
    Returns false if no more indices are available. */
int roll_ind_comb(int **ind, int n, int nmax);

/** Get the number of bins in the 1 or 2 dependent features */
void get_binx_biny(int dimensions, const struct mdt_type *mdt,
                   const char *routine, int *nbinx, int *nbiny, int *ierr);

/** Return the sum of a set. */
double get_sum(const double bin[], int nbins);

/** Return the entropy of a set. */
double entrp1(const double frq[], int nbinx);

/** Calculate the pdf p(x/independent features) summed over all
    independent features and their values except for the n_feat_fix fixed
    independent features. */
void getfrq(const struct mdt_type *mdt, const int i_feat_fix[], int n_feat_fix,
            const int i_val_fix[], int nbinx, double frq[]);

/** Return entropy of p(x/y,z,...) where y,z are the independent features.
    See pages 480-483 in Numerical Recipes for equations. */
double entrp2(double summdt, const int i_feat_fix[], const struct mdt_type *mdt,
              int n_feat_fix, int nbinx, float sumi[]);

/** Get the chi^2, etc for pdf p(x/y,z,...) */
double chisqr(double summdt, const int i_feat_fix[], const struct mdt_type *mdt,
              int n_feat_fix, int nbinx, float sumi[], double *df, double *prob,
              double *ccc, double *cramrv, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* __MDT_UTIL_H */
