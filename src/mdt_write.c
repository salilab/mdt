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
static void write_mdt_header(FILE *fp, const struct mod_mdt *mdt,
                             const struct mdt_library *mlib)
{
  int i;
  char *version = mod_short_version_get();

  fprintf(fp, "REAL OUTPUT OF PROGRAM:  ** MDT %s **\n\n\n", version);
  fprintf(fp, "Number of alignments                 : %14d\n"
          "Number of proteins in the alignments : %14d\n"
          "Number of proteins or their pairs    : %14d\n"
          "Sample size                          : %#14.6g\n\n",
          mdt->nalns, mdt->n_proteins, mdt->n_protein_pairs, mdt->sample_size);

  fprintf(fp, "FEATURES TABULATED IN THIS MULTIDIMENSIONAL TABLE:\n\n");
  fprintf(fp, "  #  FEATURE NBINS NAME\n");

  for (i = 0; i < mdt->nfeat; i++) {
    int ifeat;
    struct mod_mdt_libfeature *feat;
    ifeat = mdt->features[i].ifeat - 1;
    feat = &mlib->base.features[ifeat];
    fprintf(fp, "%3d %8d %5d %s\n", i + 1, ifeat + 1, feat->nbins, feat->name);
  }

  fprintf(fp, "\n\nSUBSET OF VALUES (WITH RESPECT TO THE BIN FILE)\n"
          "THAT ARE ACTUALLY PRESENT IN THIS FILE (last indx rolls first):\n\n"
          "  # ISTART   IEND\n");

  for (i = 0; i < mdt->nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->features[i];
    fprintf(fp, "%3d %6d %6d\n", i + 1, feat->istart, feat->iend);
  }

  fprintf(fp, "\n\nMDT TABLE START:%9d\n", mdt->nelems);

  g_free(version);
}

/** Write MDT footer to a file */
static void write_mdt_footer(FILE *fp, const struct mod_mdt *mdt)
{
  fprintf(fp, "MDT TABLE END%s\n", mdt->pdf ? ":PDF" : "");
}

/** Write MDT bin data to a file */
static void write_mdt_data(FILE *fp, const struct mod_mdt *mdt)
{
  int i;

  for (i = 0; i < mdt->nelems; i++) {
    fprintf(fp, "%#15.5g\n", mdt->bin[i]);
  }
}

/** Write out an MDT in text format. Return TRUE on success. */
gboolean mdt_write(const struct mdt *mdt, const struct mdt_library *mlib,
                   const char *filename, gboolean write_preamble, GError **err)
{
  FILE *fp;
  struct mod_file file_info;

  fp = mdt_open_file(filename, "w", &file_info, err);
  if (fp) {
    if (write_preamble) {
      write_mdt_header(fp, &mdt->base, mlib);
    }
    write_mdt_data(fp, &mdt->base);
    if (write_preamble) {
      write_mdt_footer(fp, &mdt->base);
    }
    return mdt_close_file(fp, &file_info, err);
  } else {
    return FALSE;
  }
}
