/** \file mdt_write.c   Functions to write MDTs to files.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdio.h>
#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Write MDT header to a file */
static void write_mdt_header(FILE *fp, const struct mdt_type *mdt,
                             const struct mdt_library *mlib)
{
  int i;
  char *version = get_mod_short_version();

  fprintf(fp, "REAL OUTPUT OF PROGRAM:  ** MODELLER %s **\n\n\n", version);
  fprintf(fp, "Number of alignments                 : %14d\n"
              "Number of proteins in the alignments : %14d\n"
              "Number of proteins or their pairs    : %14d\n"
              "Sample size                          : %#14.6g\n\n",
          mdt->nalns, mdt->n_proteins, mdt->n_protein_pairs, mdt->sample_size);

  fprintf(fp, "FEATURES TABULATED IN THIS MULTIDIMENSIONAL TABLE:\n\n");
  fprintf(fp, "  #  FEATURE NBINS NAME\n");

  for (i = 0; i < mdt->nfeat; i++) {
    int ifeat;
    struct mdt_feature *feat;
    ifeat = f_int1_get(&mdt->ifeat, i) - 1;
    feat = &mlib->base.features[ifeat];
    fprintf(fp, "%3d %8d %5d %s\n", i+1, ifeat+1, feat->nbins, feat->name);
  }

  fprintf(fp, "\n\nSUBSET OF VALUES (WITH RESPECT TO THE BIN FILE)\n"
          "THAT ARE ACTUALLY PRESENT IN THIS FILE (last indx rolls first):\n\n"
          "  # ISTART   IEND\n");

  for (i = 0; i < mdt->nfeat; i++) {
    fprintf(fp, "%3d %6d %6d\n", i+1, f_int1_get(&mdt->istart, i),
            f_int1_get(&mdt->iend, i));
  }

  fprintf(fp, "\n\nMDT TABLE START:%9d\n", mdt->nelems);

  free(version);
}

/** Write MDT footer to a file */
static void write_mdt_footer(FILE *fp, const struct mdt_type *mdt)
{
  fprintf(fp, "MDT TABLE END%s\n", mdt->pdf ? ":PDF" : "");
}

/** Write MDT bin data to a file */
static void write_mdt_data(FILE *fp, const struct mdt_type *mdt)
{
  int i;
  double *bin = f_double1_pt(&mdt->bin);

  for (i = 0; i < mdt->nelems; i++) {
    fprintf(fp, "%#15.5g\n", bin[i]);
  }
}

/** Write out an MDT. Return TRUE on success. */
gboolean mdt_write(const struct mdt_type *mdt, const struct mdt_library *mlib,
                   const char *filename, gboolean write_preamble, GError **err)
{
  FILE *fp;
  struct mod_file file_info;

  fp = mdt_open_file(filename, "w", &file_info, err);
  if (fp) {
    if (write_preamble) {
      write_mdt_header(fp, mdt, mlib);
    }
    write_mdt_data(fp, mdt);
    if (write_preamble) {
      write_mdt_footer(fp, mdt);
    }
    return mdt_close_file(fp, &file_info, err);
  } else {
    return FALSE;
  }
}
