/** \file mdt_section.c    Functions to handle subsections of MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <stdlib.h>
#include <math.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Given the indices of an MDT section, return the start and end indices
    into the bin array. */
static void get_mdt_section_bins(const struct mdt_type *mdt,
                                 const int indices[], int n_indices,
                                 int *istart, int *nbins, int *ierr)
{
  static const char *routine = "get_mdt_section_bins";
  int i, *indf;

  *ierr = 0;
  if (n_indices < 0 || n_indices >= mdt->nfeat) {
    modlogerror(routine, ME_VALUE, "Incorrect number of features (%d);\n"
                "must be less than the dimension of the MDT (%d)", n_indices,
                mdt->nfeat);
    *ierr = 1;
    return;
  }

  indf = mdt_start_indices(mdt);
  for (i = 0; i < n_indices; i++) {
    indf[i] = indices[i] + 1;
  }
  *istart = indmdt(indf, mdt);
  free(indmdt);
  if (*istart < 0 || *istart >= mdt->nelems) {
    modlogerror(routine, ME_INDEX, "Index %d out of range %d to %d",
                *istart, 0, mdt->nelems);
    *ierr = 1;
    return;
  }

  *nbins = 1;
  for (i = n_indices; i < mdt->nfeat; i++) {
    (*nbins) *= f_int1_get(&mdt->nbins, i);
  }
}


/** Calculate the entropy of a histogram. */
static double entropy_hist(const double x[], int n)
{
  static const float divisor = 1.0e-15, tiny = 1.0e-30;
  double area;

  area = get_sum(x, n);
  if (fabs(area) < divisor && n > 0) {
    /* the curve is probably all 0. */
      return log(n);
  } else {
    int i;
    double accum;
    accum = 0.0;
    for (i = 0; i < n; i++) {
      double probi = x[i] / area;
      if (probi > tiny) {
        accum -= probi * log(probi);
      }
    }
    return accum;
  }
}


/** Sum an MDT section. */
double mdt_section_sum(const struct mdt_type *mdt, const int indices[],
                       int n_indices, int *ierr)
{
  int istart, nbins;
  double *bin;
  *ierr = 0;
  get_mdt_section_bins(mdt, indices, n_indices, &istart, &nbins, ierr);
  if (*ierr) {
    return 0.0;
  }
  bin = f_double1_pt(&mdt->bin);
  return get_sum(&bin[istart], nbins);
}

/** Get the entropy of an MDT section. */
double mdt_section_entropy(const struct mdt_type *mdt, const int indices[],
                           int n_indices, int *ierr)
{
  int istart, nbins;
  double *bin;
  *ierr = 0;
  get_mdt_section_bins(mdt, indices, n_indices, &istart, &nbins, ierr);
  if (*ierr) {
    return 0.0;
  }
  bin = f_double1_pt(&mdt->bin);
  return entropy_hist(&bin[istart], nbins);
}
