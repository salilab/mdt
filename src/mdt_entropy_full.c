/** \file mdt_entropy_full.c  Functions to calculate entropy of an MDT.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"
#include "num_recipes.h"

/** Write an explanatory header about entropy to the log. */
static void wrhead(double hx, const struct mdt_type *mdt,
                   const struct mdt_library *mlib, double summdt, double hx0,
                   double ux0, double dfu, double chi2u, double chi2up)
{
  int i;
  char rule[80];
  for (i = 0; i < 79; i++) {
    rule[i] = '-';
  }
  rule[79] = 0;
  mod_lognote("\n\n\nEntropy output\n\n\n"
              "Features tabulated in the multidimensional table\n\n"
              "  #  FEATURE NAME");
  for (i = 0; i < mdt->nfeat; i++) {
    struct mdt_libfeature *feat;
    int ifeat = mdt->features[i].ifeat - 1;
    feat = &mlib->base.features[ifeat];
    mod_lognote("%3d   %6d %s", i + 1, ifeat + 1, feat->name);
  }

  mod_lognote("\n  The last feature is the dependent one (x).");

  mod_lognote("\n\nThe subset of values (with respect to the BIN file)"
              "\nthat are actually present in this analysis:\n\n  # NBINS");
  for (i = 0; i < mdt->nfeat; i++) {
    mod_lognote("%3d %7d", i + 1, mdt->features[i].nbins);
  }

  mod_lognote("\n\nNumber of all MDT table elements :%12d"
              "\nNumber of all MDT table points   :%12.2f", mdt->nelems,
              summdt);

  mod_logout("\n\n\nDegrees of freedom for the p(x) chi^2 test    : %14.4f\n"
             "chi^2 for comparison of p(x) with uniform pdf : %14.4f\n"
             "The significance of this chi^2                : %14.6g\n"
             "The entropy of uniform pdf, H(u)              : %12.4f\n"
             "The entropy of pdf p(x), H(x)                 : %12.4f\n"
             "U(x) = (H(u) - H(x)) / H(u)                   : %12.4f",
             dfu, chi2u, chi2up, hx0, hx, ux0);

  mod_lognote("\n\nThe significance (as measured by chi^2) and strength\n"
              "(as measured by entropy) of associations between x and\n"
              "various combinations of y,z,... :\n\n"
              "The situation with more than one independent variable is "
              "formally\n"
              "the same as if they were transformed into one super-variable "
              "with\n"
              "ny*nz different values, where ny and nz are the numbers of "
              "values\n"
              "for Y and Z variables alone.\n\n"
              "    DF       ...  degrees of freedom for the chi^2;\n"
              "    chi^2    ...  chi^2 statistic for significance of "
              "association;\n"
              "    p        ...  significance of chi^2 (small -> "
              "significant);\n"
              "    H(x/yz)  ...  weighted average entropy of the "
              "conditional pdf's\n"
              "    U(x/yz)  ...  (H(x) - H(x/yz)) / H(x),\n"
              "    u(x/yz)  ...  (H(u) - H(x/yz)) / H(u).\n");
  mod_logout("\nIndpndt feats y,z,...%9s%11s%11s%8s%8s%8s\n%s",
             "DF", "chi^2", "p", "H(x/yz)", "U(x/yz)", "u(x/yz)", rule);
}

/** Write entropy results to the log. */
static void wrres(const int i_feat_fix[], int n_feat_fix, double df,
                  double chisq, double prob, double hxy, double uxy,
                  double uuxy)
{
  int i;
  GString *str = g_string_new(NULL);
  for (i = 0; i < n_feat_fix; i++) {
    g_string_append_printf(str, "%3d", i_feat_fix[i] + 1);
  }
  g_string_truncate(str, 21);
  mod_logout("%-21s%10.2f %10.4g %10.4g %7.4f %7.4f %7.4f", str->str, df,
             chisq, prob, hxy, uxy, uuxy);
  g_string_free(str, TRUE);
}


/** Write chi^2, entropies and dependencies of pdf p(x/y,z,...) to the log.
    Return TRUE on success. */
gboolean mdt_entropy_full(const struct mdt_type *mdt,
                          const struct mdt_library *mlib, GError **err)
{
  static const double small = 1.e-8;
  const char *routine = "mdt_entropy_full";
  double summdt, dfu, expctd, chi2u, chi2up, sumfrq, hx, hx0, ux0, *frq;
  float *sumi;
  GError *tmperr = NULL;
  int nbinx, i, *i_feat_fix, n_feat_fix;

  nbinx = mdt->features[mdt->nfeat - 1].nbins;

  /* the number of points in mdt */
  summdt = get_sum(mdt->bin, mdt->nelems);

  if (summdt < small) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "%s: MDT is empty; sum over all elements = %.4g", routine,
                summdt);
    return FALSE;
  }

  /* get pdf p(x) irrespective of the values of the independent variables */
  sumi = g_malloc(sizeof(float) * nbinx);
  frq = g_malloc(sizeof(double) * nbinx);
  getfrq(mdt, NULL, 0, NULL, nbinx, frq);


  /* get its chi^2 characterization by comparison with the uniform pdf: */

  /* degrees of freedom: */
  dfu = nbinx - 1;

  /* chi^2 */
  expctd = summdt / (double)nbinx;
  chi2u = 0.;
  for (i = 0; i < nbinx; i++) {
    chi2u += (frq[i] - expctd) * (frq[i] - expctd) / expctd;
  }

  /* chi^2 significance (small is significant) */
  chi2up = gammq(0.5 * dfu, 0.5 * chi2u);

  /* get its entropy: */
  sumfrq = get_sum(frq, nbinx);
  sumfrq = (sumfrq < small ? small : sumfrq);
  for (i = 0; i < nbinx; i++) {
    frq[i] /= sumfrq;
  }
  hx = entrp1(frq, nbinx);
  if (hx < small) {
    mod_logwarning(routine, "Entropy too small for division; changed to %e",
                   small);
    hx = small;
  }

  /* entropy of the equivalent uniform distribution */
  hx0 = -log(1. / (double)nbinx);
  ux0 = (hx0 - hx) / hx0;

  /* write out the header: */
  wrhead(hx, mdt, mlib, summdt, hx0, ux0, dfu, chi2u, chi2up);

  /* for each number of independent variables from 1 to nfeat-1 */
  for (n_feat_fix = 1; n_feat_fix < mdt->nfeat; n_feat_fix++) {
    int ncomb2;

    /* the number of different combinations of n_feat_fix features from the
       set of numb_features-1 possible features */
    ncomb2 = nperm(mdt->nfeat - 1) / (nperm(n_feat_fix)
                                      * nperm(mdt->nfeat - 1 - n_feat_fix));

    /* generate all combinations of independent features
       (there are n_feat_fix of these features, ncomb2 combinations); */
    i_feat_fix = NULL;
    while (roll_ind_comb(&i_feat_fix, n_feat_fix, mdt->nfeat - 1)) {
      double hxy, chisq, uxy, uuxy, df, prob, ccc, cramrv;

      /* get the conditional entropy of pdf p(x/y,z,...): */
      hxy = entrp2(summdt, i_feat_fix, mdt, n_feat_fix, nbinx, sumi);

      /* get the chi^2, etc for pdf p(x/y,z,...): */
      chisq = chisqr(summdt, i_feat_fix, mdt, n_feat_fix, nbinx, sumi, &df,
                     &prob, &ccc, &cramrv, &tmperr);
      if (tmperr) {
        break;
      }

      /* get the uncertainty coefficient of pdf p(x/y,z,...) (comparison
         of p(x/yz) with p(x)): */
      uxy = (hx - hxy) / hx;
      /* get the equivalent measure for comparison with uniform pdf: */
      uuxy = (hx0 - hxy) / hx0;

      /* write out: */
      wrres(i_feat_fix, n_feat_fix, df, chisq, prob, hxy, uxy, uuxy);
    }
    g_free(i_feat_fix);
    if (tmperr) {
      break;
    }
  }
  g_free(frq);
  g_free(sumi);
  if (tmperr) {
    g_propagate_error(err, tmperr);
    return FALSE;
  } else {
    return TRUE;
  }
}
