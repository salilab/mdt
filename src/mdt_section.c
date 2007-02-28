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
  free(indf);
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

/** Returns arg1/arg2 unless arg2 is so small that division
    cannot be done accurately. In that case, 0 is returned. */
static float divide2(float arg1, float arg2)
{
  static const float divisor = 1.0e-15;

  if (fabs(arg2) < divisor) {
    return 0.0;
  } else {
    return arg1 / arg2;
  }
}


/** Returns the difference (in degrees) of two angles
    (in degrees): ANG2-ANG1.  The difference is in the interval
    from -180 -- 180 degrees. It is defined as the shortest path
    around the circle from ANG1 to ANG2 (clockwise is +). */
static float diffdeg(float ang2, float ang1)
{
  static const float pi2degr = 360.0, pidegr = 180.0;
  float a1, a2, d;

  a1 = fmod(ang1, pi2degr);
  a2 = fmod(ang2, pi2degr);
  d = a2 - a1;
  if (d < (-pidegr)) {
    d += pi2degr;
  }
  if (d > pidegr) {
    d -= pi2degr;
  }
  return d;
}


/** Return the average and standard deviation of the distribution
    defined by f(i). Works only for aperiodic features. */
static void hist_avrstdev(const double f[], int n, double x0, double dx,
                          double *avr, double *std)
{
  float a, s, x, norm;
  int i;

  a = norm = 0.;
  x = x0;
  for (i = 0; i < n; i++) {
    norm += f[i];
    a += f[i] * x;
    x += dx;
  }

  *avr = divide2(a, norm);

  /* although it could do both in the same loop
     (cf, algebra_routines.F:avrstdev), it is probably more
     numerically stable to do it this way */
  s = 0.;
  x = x0;
  for (i = 0; i < n; i++) {
    s += f[i] * pow(x - *avr, 2);
    x += dx;
  }

  /* standard deviation: */
  *std = sqrt(divide2(s, norm));
}


/** Return the average and standard deviation of the distribution
    defined by f(i). Works only for 360^o periodic features between
    -180 and 180^o (degrees only, no radians) that have single
    relatively sharp peaks (such as bond angles and improper dihedral
    angles in proteins). you can check the results from here against
    the analytical model fits obtained with ASGL. Bond angles
    (0 .. 180^o) don't need complication, but it does not hurt. */
static void hist_avrstdev_deg(const double f[], int n, float x0, float dx,
                              double *avr, double *std)
{
  static const float pi2degr = 360.0, pidegr = 180.0;
  float a1, a2, s, x, norm1, norm2;
  int i;

  /* need to employ a trick: average positive values, average negative
     values then average the two, taking into account periodicity */
  a1 = a2 = norm1 = norm2 = 0.0;
  x = x0;
  for (i = 0; i < n; i++) {
    if (x <= 0.) {
      a1 += f[i] * x;
      norm1 += f[i];
    } else {
      a2 += f[i] * x;
      norm2 += f[i];
    }
    x += dx;
  }

  a1 = divide2(a1, norm1);
  a2 = divide2(a2, norm2);

  /* are a1 and a2 closer if they are "connected" via 0 or 180 degrees? */
  if (a2 - a1 > pidegr) {
    a1 += pi2degr;
  }

  *avr = divide2(norm1 * a1 + norm2 * a2, norm1 + norm2);

  /* although it could do both in the same loop
     (cf, algebra_routines.F:avrstdev), is is probably more
     numerically stable to do it this way */
  s = 0.;
  x = x0;
  for (i = 0; i < n; i++) {
    s += f[i] * pow(diffdeg(x, *avr), 2);
    x += dx;
  }

  /* standard deviation: */
  *std = sqrt(divide2(s, norm1 + norm2));
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

/** Get the mean and standard deviation of an MDT section. */
void mdt_section_meanstdev(const struct mdt_type *mdt,
                           const struct mdt_library *mlib, const int indices[],
                           int n_indices, double *mean, double *stdev,
                           int *ierr)
{
  int istart, nbins, ifeat;
  mbool periodic;
  double *bin, dx, x0;

  *ierr = 0;
  get_mdt_section_bins(mdt, indices, n_indices, &istart, &nbins, ierr);
  if (*ierr) {
    return;
  }
  bin = f_double1_pt(&mdt->bin);
  ifeat = f_int1_get(&mdt->ifeat, mdt->nfeat - 1);
  periodic = mdt_feature_is_periodic(ifeat);

  /* histogram bin size in real units: */
  dx = f_float2_get(&mlib->rang2, 0, ifeat-1)
       - f_float2_get(&mlib->rang1, 0, ifeat-1);
  /* position of the center of the first bin in real units: */
  x0 = f_float2_get(&mlib->rang1, 0, ifeat-1) + 0.5 * dx;

  if (periodic) {
    hist_avrstdev_deg(&bin[istart], nbins, x0, dx, mean, stdev);
  } else {
    hist_avrstdev(&bin[istart], nbins, x0, dx, mean, stdev);
  }
}