/** \file num_recipes.c    Numerical functions, largely from Numerical Recipes.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
 */

#include <assert.h>
#include <math.h>

#include "modeller.h"
#include "num_recipes.h"

/** Return the number of combinations of n. */
int nperm(int n)
{
  int np = 1, i;
  for (i = 2; i <= n; i++) {
    np *= i;
  }
  return np;
}

/** Return the logarithm of the gamma function. */
static float gammln(float xx)
{
  float x, tmp, ser;
  static const float cof[6] = { 76.18009173, -86.50532033, 24.01409822,
    -1.231739516, 0.120858003e-2, -0.536382e-5
  };
  static const float stp = 2.50662827465;
  int j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for (j = 0; j < 6; j++) {
    x += 1.0;
    ser += cof[j] / x;
  }
  return -tmp + log(stp * ser);
}

/** Generate a continued fraction used by gammq(). */
static void gcf(float *gammcf, float a, float x, float *gln)
{
  static const char *routine = "gcf";
  static const int itmax = 100;
  static const float eps = 3.0e-7;
  int n;
  float gold = 0.0, g, fac = 1.0, b1 = 1.0, b0 = 0.0, anf, ana, an, a1,
      a0 = 1.0;

  *gln = gammln(a);
  a1 = x;
  for (n = 1; n <= itmax; n++) {
    an = (float)n;
    ana = an - a;
    a0 = (a1 + a0 * ana) * fac;
    b0 = (b1 + b0 * ana) * fac;
    anf = an * fac;
    a1 = x * a0 + anf * a1;
    b1 = x * b0 + anf * b1;
    if (a1) {
      fac = 1.0 / a1;
      g = b1 * fac;
      if (fabs((g - gold) / g) < eps) {
        *gammcf = exp(-x + a * log(x) - (*gln)) * g;
        return;
      }
      gold = g;
    }
  }
  mod_logwarning(routine, "a too large or itmax too small");
  *gammcf = 0.;
}

/** Generate a series used by gammq(). */
static void gser(float *gamser, float a, float x, float *gln)
{
  static const char *routine = "gser";
  static const int itmax = 100;
  static const float eps = 3.0e-7;

  *gln = gammln(a);
  assert(x >= 0.);
  if (x == 0.) {
    *gamser = 0.;
    return;
  } else {
    int n;
    float sum, del, ap;

    ap = a;
    del = sum = 1. / a;
    for (n = 1; n <= itmax; n++) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (fabs(del) < fabs(sum) * eps) {
        *gamser = sum * exp(-x + a * log(x) - (*gln));
        return;
      }
    }
    mod_logwarning(routine, "a too large or itmax too small");
    *gamser = 0.;
  }
}


/** Return the complement of incomplete gamma function. */
float gammq(float a, float x)
{
  float gln;
  assert(x >= 0. && a > 0.);

  if (x < a + 1.) {
    float gamser;
    gser(&gamser, a, x, &gln);
    return 1.0 - gamser;
  } else {
    float gammcf;
    gcf(&gammcf, a, x, &gln);
    return gammcf;
  }
}
