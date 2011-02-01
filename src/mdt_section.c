/** \file mdt_section.c    Functions to handle subsections of MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <stdlib.h>
#include <math.h>
#include "modeller.h"
#include "mdt.h"
#include "mdt_feature.h"
#include "util.h"

/** Given the indices of an MDT section, return the start and end indices
    into the bin array. */
static gboolean get_mdt_section_bins(const struct mod_mdt *mdt,
                                     const int indices[], int n_indices,
                                     int *istart, int *nbins, GError **err)
{
  static const char *routine = "get_mdt_section_bins";
  int i, *indf;

  if (n_indices < 0 || n_indices >= mdt->nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: Incorrect number of features (%d);\n"
                "must be less than the dimension of the MDT (%d)", routine,
                n_indices, mdt->nfeat);
    return FALSE;
  }

  indf = mdt_start_indices(mdt);
  for (i = 0; i < n_indices; i++) {
    indf[i] = indices[i] + 1;
  }
  *istart = indmdt(indf, mdt);
  free(indf);
  if (*istart < 0 || *istart >= mdt->nelems) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                "%s: Index %d out of range %d to %d", routine, *istart, 0,
                mdt->nelems);
    return FALSE;
  }

  *nbins = 1;
  for (i = n_indices; i < mdt->nfeat; i++) {
    (*nbins) *= mdt->features[i].nbins;
  }
  return TRUE;
}


/** Calculate the entropy of a histogram. */
static double entropy_hist(const struct mod_mdt *mdt, int offset, int n)
{
  static const float divisor = 1.0e-15, tiny = 1.0e-30;
  int i;
  double area = 0.0;

  for (i = 0; i < n; ++i) {
    area += mod_mdt_bin_get(mdt, offset + i);
  }
  if (fabs(area) < divisor && n > 0) {
    /* the curve is probably all 0. */
    return log(n);
  } else {
    int i;
    double accum;
    accum = 0.0;
    for (i = 0; i < n; i++) {
      double probi = mod_mdt_bin_get(mdt, offset + i) / area;
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
static void hist_avrstdev(const struct mod_mdt *mdt, int offset, int n,
                          double x0, double dx, double *avr, double *std)
{
  float a, s, x, norm;
  int i;

  a = norm = 0.;
  x = x0;
  for (i = 0; i < n; i++) {
    double binval = mod_mdt_bin_get(mdt, offset + i);
    norm += binval;
    a += binval * x;
    x += dx;
  }

  *avr = divide2(a, norm);

  /* although it could do both in the same loop
     (cf, algebra_routines.F:avrstdev), it is probably more
     numerically stable to do it this way */
  s = 0.;
  x = x0;
  for (i = 0; i < n; i++) {
    double binval = mod_mdt_bin_get(mdt, offset + i);
    s += binval * pow(x - *avr, 2);
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
static void hist_avrstdev_deg(const struct mod_mdt *mdt, int offset, int n,
                              float x0, float dx, double *avr, double *std)
{
  static const float pi2degr = 360.0, pidegr = 180.0;
  float a1, a2, s, x, norm1, norm2;
  int i;

  /* need to employ a trick: average positive values, average negative
     values then average the two, taking into account periodicity */
  a1 = a2 = norm1 = norm2 = 0.0;
  x = x0;
  for (i = 0; i < n; i++) {
    double binval = mod_mdt_bin_get(mdt, offset + i);
    if (x <= 0.) {
      a1 += binval * x;
      norm1 += binval;
    } else {
      a2 += binval * x;
      norm2 += binval;
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
    double binval = mod_mdt_bin_get(mdt, offset + i);
    s += binval * pow(diffdeg(x, *avr), 2);
    x += dx;
  }

  /* standard deviation: */
  *std = sqrt(divide2(s, norm1 + norm2));
}


/** Sum an MDT section. */
double mdt_section_sum(const struct mod_mdt *mdt, const int indices[],
                       int n_indices, GError **err)
{
  int istart, nbins, i;
  double sum = 0.0;
  if (!get_mdt_section_bins(mdt, indices, n_indices, &istart, &nbins, err)) {
    return 0.0;
  }
  for (i = 0; i < nbins; ++i) {
    sum += mod_mdt_bin_get(mdt, istart + i);
  }
  return sum;
}

/** Get the entropy of an MDT section. */
double mdt_section_entropy(const struct mod_mdt *mdt, const int indices[],
                           int n_indices, GError **err)
{
  int istart, nbins;
  if (!get_mdt_section_bins(mdt, indices, n_indices, &istart, &nbins, err)) {
    return 0.0;
  }
  return entropy_hist(mdt, istart, nbins);
}

/** Get the mean and standard deviation of an MDT section. */
void mdt_section_meanstdev(const struct mod_mdt *mdt,
                           const struct mdt_library *mlib, const int indices[],
                           int n_indices, double *mean, double *stdev,
                           GError **err)
{
  int istart, nbins, ifeat;
  gboolean periodic;
  struct mod_mdt_libfeature *feat;
  double dx, x0;

  if (!get_mdt_section_bins(mdt, indices, n_indices, &istart, &nbins, err)) {
    *mean = *stdev = 0.;
    return;
  }
  ifeat = mdt->features[mdt->nfeat - 1].ifeat;
  feat = &mlib->base.features[ifeat - 1];
  periodic = mdt_feature_periodic_get(mlib, ifeat);

  /* histogram bin size in real units: */
  dx = feat->bins[0].rang2 - feat->bins[0].rang1;
  /* position of the center of the first bin in real units: */
  x0 = feat->bins[0].rang1 + 0.5 * dx;

  if (periodic) {
    hist_avrstdev_deg(mdt, istart, nbins, x0, dx, mean, stdev);
  } else {
    hist_avrstdev(mdt, istart, nbins, x0, dx, mean, stdev);
  }
}
