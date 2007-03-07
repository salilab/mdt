/** \file mdt_write_asgl.c    Functions to write input files for ASGL.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Write the raw data to be plotted with ASGL. */
static void wrdata(const char *datfil, int dimensions, const double bin[],
                   int nbinx, int nbiny, int *ierr)
{
  FILE *fp;
  struct mod_file file_info;
  fp = open_file(datfil, "w", &file_info);
  if (fp) {
    int i;
    for (i = 0; i < nbinx * nbiny; i++) {
      fprintf(fp, "%#14.5g\n", bin[i]);
    }
    close_file(fp, &file_info, ierr);
  } else {
    *ierr = 1;
  }
}

/** Get the symbol name for a bin. */
static char *get_mdt_symb(const struct mdt_type *mdt,
                          const struct mdt_library *mlib,
                          int ifeat, int ibin, int ndecimal)
{
  /* For type 3, generate symbol from range data */
  if (mlib->itsymb[ifeat] == 3) {
    if (ibin == mlib->ndimen[ifeat] - 1) {
      return g_strdup("U");
    } else {
      float rang1 = f_float2_get(&mlib->rang1, ibin, ifeat);
      if (ndecimal > 0) {
        char *fmt = g_strdup_printf("%%.%df", ndecimal);
        char *str = g_strdup_printf(fmt, rang1);
        g_free(fmt);
        return str;
      } else {
        return g_strdup_printf("%d", (int)rang1);
      }
    }
  } else {
    return mdt_symb_get(mdt, mlib, ifeat, ibin);
  }
}

/** Output the ASGL commands to make a single plot. */
static void appasgl(FILE *fp, const struct mdt_type *mdt,
                    const struct mdt_library *mlib, const char *datfil,
                    int ipos, const int indf[], int dimensions,
                    int every_x_numbered, int every_y_numbered, int nbinx,
                    int nbiny, const char *text, const char *plot_type,
                    int x_decimal, int y_decimal, double sum)
{
  char *featnam;
  int i, ifeat, *ifeatpt, itsymbx, itsymby;
  ifeatpt = f_int1_pt(&mdt->ifeat);

  ifeat = ifeatpt[mdt->nfeat - 1];
  itsymbx = mlib->itsymb[ifeat - 1];
  ifeat = ifeatpt[mdt->nfeat > 1 ? mdt->nfeat - 2 : 0];
  itsymby = mlib->itsymb[ifeat - 1];

  fputs("# -------------------------------------------------\n", fp);
  if (dimensions == 1) {
    fprintf(fp, "READ_TABLE FILE '%s'\n", datfil);
  } else {
    fprintf(fp, "READ_DPLOT DPLOT_ORIENTATION = 'YX', FILE = '%s'\n", datfil);
  }

  /* possibly make a smart correction to explicitly specified nevryx
     nevryy: only if it is 1 and the automatic ranges are used: */
  if (every_x_numbered == 1 && itsymbx == 3) {
    every_x_numbered = MAX(1, nbinx / 8 + 1);
  }
  if (every_y_numbered == 1 && itsymby == 3) {
    every_y_numbered = MAX(1, nbiny / 8 + 1);
  }

  if (dimensions == 1) {
    fprintf(fp, "SET X_LABEL_STYLE = 1\n"
                "SET X_TICK = 1 1 -999\n"
                "SET POSITION %d 0\n"
                "SET XY_COLUMNS 0 1\n", ipos);
    if (itsymbx == 2 || itsymbx == 3) {
      fprintf(fp, "SET WORLD_WINDOW 0.5 0 %7.1f %7.1f\n", nbinx + 1.5, -999.);
    } else {
      fprintf(fp, "SET WORLD_WINDOW 0 0 %5d %5d\n", nbinx + 1, -999);
    }
  } else {
    fprintf(fp, "SET X_LABEL_STYLE = 1\n"
                "SET POSITION %5d 1\n"
                "SET Y_LABEL_STYLE = 1\n"
                "SET Y_TICK_LABEL = 1\n", ipos);
    if (itsymbx == 2 || itsymbx == 3) {
      fputs("SET Y_TICK = 0.5 1 -999\nSET X_TICK = 0.5 1 -999\n", fp);
    } else {
      fprintf(fp, "SET Y_TICK = 1 1 %3d\nSET X_TICK = 1 1 %3d\n", nbiny, nbinx);
    }
    fprintf(fp, "SET WORLD_WINDOW 0 0 %5d %5d\n", nbinx + 1, nbiny + 1);
    ifeat = ifeatpt[mdt->nfeat - 2];
    fputs("SET Y_LABELS", fp);
    for (i = 0; i < nbiny; i += every_y_numbered) {
      char *symb = get_mdt_symb(mdt, mlib, ifeat - 1, i, y_decimal);
      fprintf(fp, " '%s'", symb);
      g_free(symb);
    }
    fputs("\n", fp);
  }

  /* BAR_XSHIFT does not have any effect in DPLOT, so this part is the same
     for both HIST2D and DPLOT; the first tick must be labelled: */
  if (itsymbx == 2 ||  itsymbx == 3) {
    fprintf(fp, "SET X_TICK_LABEL = 1 %3d, BAR_XSHIFT = 0.5\n",
            every_x_numbered);
  } else {
    fprintf(fp, "SET X_TICK_LABEL = 1 %3d\n", every_x_numbered);
  }

  fputs(text, fp);
  fputs("\n", fp);

  ifeat = ifeatpt[mdt->nfeat - 1];
  fputs("SET X_LABELS", fp);
  for (i = 0; i < nbinx; i += every_x_numbered) {
    char *symb = get_mdt_symb(mdt, mlib, ifeat - 1, i, x_decimal);
    fprintf(fp, " '%s'", symb);
    g_free(symb);
  }
  fputs("\nWORLD\nAXES2D\nRESET_CAPTIONS\n", fp);

  fprintf(fp, "CAPTION CAPTION_POSITION 1, ;\n"
              "     CAPTION_TEXT '%.1f POINTS'\n", sum);

  ifeat = ifeatpt[mdt->nfeat - 1];
  featnam = mdt_library_featnam_get(mlib, ifeat - 1);
  fprintf(fp, "CAPTION CAPTION_POSITION 2, ;\n"
              "     CAPTION_TEXT '%s'\n", featnam);
  g_free(featnam);

  if (dimensions == 1) {
    fputs("CAPTION CAPTION_POSITION 3, ;\n"
          "     CAPTION_TEXT 'FREQUENCY'\n", fp);
  } else {
    ifeat = ifeatpt[mdt->nfeat - 2];
    featnam = mdt_library_featnam_get(mlib, ifeat - 1);
    fprintf(fp, "CAPTION CAPTION_POSITION 3, ;\n"
                "     CAPTION_TEXT '%s'\n", featnam);
    g_free(featnam);
  }

  for (i = mdt->nfeat - dimensions - 1; i >= 0; i--) {
    char *symb;
    ifeat = ifeatpt[i];
    featnam = mdt_library_featnam_get(mlib, ifeat - 1);
    symb = mdt_symb_get(mdt, mlib, ifeat - 1, indf[i] - 1);
    fprintf(fp, "CAPTION CAPTION_POSITION 1, ;\n"
                "     CAPTION_TEXT '%s : %s'\n", strlen(symb) == 0 ? "u" : symb,
                featnam);
    g_free(symb);
    g_free(featnam);
  }

  if (dimensions == 1) {
    fputs(strcmp(plot_type, "PLOT2D") == 0 ? "PLOT2D\n" : "HIST2D\n", fp);
  } else {
    fputs("DPLOT\n", fp);
  }
}

/** Write the ASGL .top file */
static void write_script_file(FILE *fp, const struct mdt_type *mdt,
                              const struct mdt_library *mlib,
                              int dimensions, int nbinx, int nbiny,
                              double plot_density_cutoff, const char *asglroot,
                              int plots_per_page, int plot_position,
                              int every_x_numbered, int every_y_numbered,
                              const char *text, const char *plot_type,
                              int x_decimal, int y_decimal, int *ierr)
{
  int *indf, nbins, npage, nhist, nempty, iposc, idrawn;
  double *bin;

  bin = f_double1_pt(&mdt->bin);
  indf = mdt_start_indices(mdt);
  nbins = nbinx * nbiny;
  /* initialize plotting counters */
  npage = 1;
  nhist = nempty = iposc = idrawn = 0;

  do {
    int i1;
    double sum;
    /* one more attempted plot */
    nhist++;
    i1 = indmdt(indf, mdt);

    /* is it dense enough to be worth processing at all? */
    sum = get_sum(&bin[i1], nbins);
    if (sum >= plot_density_cutoff) {
      int ipos;
      char *datfil;
      idrawn++;

      /* position of the plot in the ASGL convention: */
      ipos = (iposc % plots_per_page) + plot_position;
      iposc++;

      /* get the file name for the numbers: */
      datfil = g_strdup_printf("%s.%d", asglroot, nhist);

      /* append TOP commands to the ASGL file: */
      appasgl(fp, mdt, mlib, datfil, ipos, indf, dimensions, every_x_numbered,
              every_y_numbered, nbinx, nbiny, text, plot_type, x_decimal,
              y_decimal, sum);

      /* write the numbers file */
      wrdata(datfil, dimensions, &bin[i1], nbinx, nbiny, ierr);
      g_free(datfil);

      /* new page?
         substitute nhist with ipos if you want the page jump to be decided
         by the actually drawn histograms not by the number of attempted
         plots which include those with too little points to be actually
         drawn. */
      if (ipos == (plot_position + plots_per_page - 1)) {
        fputs("NEW_PAGE\n", fp);
        iposc = 0;
        npage++;
      }
    } else {
      nempty++;
    }
  } while (roll_ind(indf, f_int1_pt(&mdt->istart), f_int1_pt(&mdt->iend),
                    mdt->nfeat - dimensions) && *ierr == 0);

  g_free(indf);
  if (*ierr) {
    return;
  }

  /* final counter processing: */
  if (idrawn == 0) {
    npage = 0;
  }

  /* Finish with a nice message: */
  modlognote("Number of attempted plots : %d\n"
             "Number of drawn plots     : %d\n"
             "Number of empty plots     : %d\n"
             "Number of pages           : %d", nhist, idrawn, nempty, npage);
}

/** Write input files to plot the given MDT with ASGL. */
void mdt_write_asgl(const struct mdt_type *mdt, const struct mdt_library *mlib,
                    const char *asglroot, const char *text, int dimensions,
                    int every_x_numbered, int every_y_numbered,
                    double plot_density_cutoff, int plots_per_page,
                    int plot_position, const char *plot_type, int x_decimal,
                    int y_decimal, int *ierr)
{
  const static char *routine = "mdt_write_asgl";
  char *topfile;
  int nbinx, nbiny;
  FILE *fp;
  struct mod_file file_info;

  *ierr = 0;
  get_binx_biny(dimensions, mdt, routine, &nbinx, &nbiny, ierr);
  if (*ierr != 0) {
    return;
  }

  topfile = g_strdup_printf("%s.top", asglroot);
  fp = open_file(topfile, "w", &file_info);
  g_free(topfile);

  if (fp) {
    write_script_file(fp, mdt, mlib, dimensions, nbinx, nbiny,
                      plot_density_cutoff, asglroot, plots_per_page,
                      plot_position, every_x_numbered, every_y_numbered, text,
                      plot_type, x_decimal, y_decimal, ierr);
    close_file(fp, &file_info, ierr);
  } else {
    *ierr = 1;
  }
}
