/** \file util.h           MDT utility functions.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#ifndef __MDT_UTIL_H
#define __MDT_UTIL_H

#include <glib.h>
#include <stdio.h>
#include "mdt_config.h"
#include "modeller.h"
#include "mod_hdf5.h"
#include "mdt_error.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
MDTDLLLOCAL
int *mdt_start_indices(const struct mod_mdt *mdt);

/** Calculate the weights in the smoothing procedure for combining the
    a priori pdf with the experimental pdf. */
MDTDLLLOCAL
void weights(float weight, int nbins, float norm, float *w1, float *w2);

/** Return the index in the MDT pdf of the point whose feature bin
    indices are indf. */
MDTDLLLOCAL
int indmdt(const int *indf, const struct mod_mdt *mdt);

/** Return the index in the MDT pdf of the point whose feature bin
    indices are indf, using the stride and istart arrays. */
MDTDLLLOCAL
int indmdt_full(const int *indf, const int stride[], int nfeat,
                const int istart[]);

/** Update the indices for the next point in the MDT. Return false if no
    more points are available. */
MDTDLLLOCAL
int roll_ind(int indf[], const int istart[], const int iend[], int nfeat);

/** Update the indices for the next point in the MDT. Return false if no
    more points are available. */
MDTDLLLOCAL
int roll_ind_mdt(int indf[], const struct mod_mdt *mdt, int nfeat);

/** Like roll_ind(), but only for the selected inds[n_inds] features */
MDTDLLLOCAL
int roll_inds(int indf[], const struct mod_mdt *mdt, const int inds[],
              int n_inds);

/** Roll n indices in ind so that all combinations of n different indices
    are generated, where the possible values for each index are 1+x to nmax-y.

    For example, for n=3:

    for (i = 0, i < nmax-2; i++) {
      for (j = i+1, j < nmax-1; j++) {
        for (k = j+1, k < nmax; k++) {

    *ind should be NULL on the first call to this routine, and it will be
    initialized (the user should free it when finished).
    Returns false if no more indices are available. */
MDTDLLLOCAL
int roll_ind_comb(int **ind, int n, int nmax);

/** Get the number of bins in the 1 or 2 dependent features. Return TRUE on
    success. */
MDTDLLLOCAL
gboolean get_binx_biny(int dimensions, const struct mod_mdt *mdt,
                       const char *routine, int *nbinx, int *nbiny,
                       GError **err);

/** Return the sum of a set. */
MDTDLLLOCAL
double get_sum(const double bin[], int nbins);

/** Return the sum of an MDT. */
MDTDLLLOCAL
double get_mdt_sum(const struct mod_mdt *mdt);

/** Return the entropy of a set. */
MDTDLLLOCAL
double entrp1(const double frq[], int nbinx);

/** Calculate the pdf p(x/independent features) summed over all
    independent features and their values except for the n_feat_fix fixed
    independent features. */
MDTDLLLOCAL
void getfrq(const struct mod_mdt *mdt, const int i_feat_fix[], int n_feat_fix,
            const int i_val_fix[], int dimensions, int nbinx, double frq[]);

/** Return entropy of p(x/y,z,...) where y,z are the independent features.
    See pages 480-483 in Numerical Recipes for equations. */
MDTDLLLOCAL
double entrp2(double summdt, const int i_feat_fix[], const struct mod_mdt *mdt,
              int n_feat_fix, int nbinx, float sumi[]);

/** Get the chi^2, etc for pdf p(x/y,z,...) */
MDTDLLLOCAL
double chisqr(double summdt, const int i_feat_fix[], const struct mod_mdt *mdt,
              int n_feat_fix, int nbinx, float sumi[], double *df, double *prob,
              double *ccc, double *cramrv, GError **err);

/** Make the stride array for faster indmdt lookup */
MDTDLLLOCAL
void make_mdt_stride(struct mod_mdt *mdt);

/** Make the stride array from the nbins array, and return the size of
    the MDT. */
MDTDLLLOCAL
int make_mdt_stride_full(const int nbins[], int nfeat, int stride[]);

/** Open a file, and uncompress it if necessary. */
MDTDLLLOCAL
struct mod_file *mdt_open_file(const char *path, const char *mode,
                               GError **err);

/** Close an open file, and do any other necessary tidy-up if it was
    compressed. The initial value of err is used (if an error was already
    set, it is not modified, but emergency cleanup is done here). */
MDTDLLLOCAL
gboolean mdt_close_file(struct mod_file *fh, GError **err);

/** Convert a Modeller error code into a GError */
MDTDLLLOCAL
void handle_modeller_error(GError **err);

/** Check that a given feature type is within range. Return FALSE and set the
    error indicator if it is not. */
MDTDLLLOCAL
gboolean check_feature_type(int ifeat, const struct mdt_library *mlib,
                            GError **err);

/** Get the position in the MDT bin array of the given set of indices. Return
    TRUE on success. */
MDTDLLLOCAL
gboolean get_bin_index(const struct mod_mdt *mdt, const int indices[],
                       int n_indices, int *bin_index, GError **err);

/** Do some basic setup of an MDT's features. Return TRUE on success. */
MDTDLLLOCAL
gboolean mdt_setup(struct mdt *mdt, const struct mdt_library *mlib,
                   GError **err);

/** Get the HDF5 datatype for this MDT. */
MDTDLLLOCAL
hid_t mdt_get_hdf5_type(const struct mod_mdt *mdt);

#ifdef WIN32
/** Our local implementation of the error function */
MDTDLLLOCAL
double erf(double x);

/** Our local implementation of the complementary error function */
MDTDLLLOCAL
double erfc(double x);
#endif

G_END_DECLS

#endif  /* __MDT_UTIL_H */
