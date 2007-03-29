/** \file mdt_normalize.c  Functions to normalize MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get the ranges of the last one or two features */
static float get_dxdy(const float dx_dy[], const struct mdt_type *mdtin,
                      int dimensions, const struct mdt_library *mlib)
{
  static const float undefined = -999;
  float dx;

  if (dx_dy[0] == undefined) {
    struct mdt_feature *feat;
    int ifeat = f_int1_get(&mdtin->ifeat, mdtin->nfeat-1) - 1;
    feat = &mlib->base.features[ifeat];
    dx = feat->bins[0].rang2 - feat->bins[0].rang1;
  } else {
    dx = dx_dy[0];
  }

  if (dimensions == 2) {
    float dy;
    if (dx_dy[1] == undefined) {
      struct mdt_feature *feat;
      int ifeat = f_int1_get(&mdtin->ifeat, mdtin->nfeat-2) - 1;
      feat = &mlib->base.features[ifeat];
      dy = feat->bins[0].rang2 - feat->bins[0].rang1;
    } else {
      dy = dx_dy[1];
    }
    return dx * dy;
  } else {
    return dx;
  }
}

/** Normalize over nbins bins */
static void do_normalize(const struct mdt_type *mdtin, struct mdt_type *mdtout,
                         int indf[], float dxdy, gboolean to_zero, int nbins,
                         int nfeat)
{
  static const float divisor = 1e-15;
  double *in_bin = f_double1_pt(&mdtin->bin);
  double *out_bin = f_double1_pt(&mdtout->bin);
  int *istart = f_int1_pt(&mdtin->istart);
  int *iend = f_int1_pt(&mdtin->iend);

  do {
    double norm;
    int i1 = indmdt(indf, mdtin);
    int i2 = i1 + nbins;

    norm = get_sum(in_bin, nbins) * dxdy;
    if (norm > divisor) {
      int i;
      for (i = i1; i < i2; i++) {
        out_bin[i] = in_bin[i] / norm;
      }
    } else if (to_zero) {
      int i;
      for (i = i1; i < i2; i++) {
        out_bin[i] = 0.0;
      }
    } else {
      int i;
      for (i = i1; i < i2; i++) {
        out_bin[i] = 1.0 / ((float)nbins * dxdy);
      }
    }

  /* roll the indices of the "constant" features one forward: */
  } while (roll_ind(indf, istart, iend, nfeat));
}


/** Normalize an MDT. Return TRUE on success. */
gboolean mdt_normalize(const struct mdt_type *mdtin, struct mdt_type *mdtout,
                       const struct mdt_library *mlib, int dimensions,
                       const float dx_dy[], int n_dx_dy, gboolean to_zero,
                       gboolean to_pdf, GError **err)
{
  static const char *routine = "mdt_normalize";
  float dxdy;
  int nbins, nbinx, nbiny, *indf;

  if (n_dx_dy != dimensions) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: dx_dy must contain %d elements, "
                "to agree with 'dimensions'.", routine, dimensions);
    return FALSE;
  }
  if (!get_binx_biny(dimensions, mdtin, routine, &nbinx, &nbiny, err)) {
    return FALSE;
  }
  nbins = nbinx * nbiny;
  if (to_pdf) {
    dxdy = get_dxdy(dx_dy, mdtin, dimensions, mlib);
  } else {
    dxdy = 1.0;
  }

  copy_mdt(mdtin, mdtout);

  modlognote("%s______> to_pdf        : %d", routine, to_pdf);
  modlognote("%s______> dimensions    : %d", routine, dimensions);
  modlognote("%s______> dx*dy         : %10.4f", routine, dxdy);
  modlognote("%s______> to_zero       : %d", routine, to_zero);

  indf = mdt_start_indices(mdtin);
  do_normalize(mdtin, mdtout, indf, dxdy, to_zero, nbins,
               mdtin->nfeat - dimensions);
  free(indf);
  mdtout->pdf = 1;
  return TRUE;
}
