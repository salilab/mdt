/** \file mdt_read.c       Functions to read MDTs from text files.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include <stdio.h>
#include <string.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Read a single line of text from the file. Return TRUE on success. */
static gboolean mdt_file_read_line(FILE *fp, GString *str, int *eof,
                                   GError **err)
{
  int ierr;
  mod_file_read_line(fp, str, &ierr, eof);
  if (ierr != 0) {
    handle_modeller_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Read an integer value from a line of text. Return TRUE on success. */
static gboolean read_int_line(const char *text, const char *starttext,
                              int intoffset, int *output, GError **err)
{
  if (strncmp(text, starttext, strlen(starttext)) == 0) {
    if (strlen(text) <= intoffset
        || sscanf(&text[intoffset], "%14d", output) != 1) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                  "Was expecting to read an integer from file line: %s", text);
      return FALSE;
    }
  }
  return TRUE;
}

/** Read a double value from a line of text. Return TRUE on success. */
static gboolean read_double_line(const char *text, const char *starttext,
                                 int intoffset, double *output, GError **err)
{
  if (strncmp(text, starttext, strlen(starttext)) == 0) {
    if (strlen(text) <= intoffset
        || sscanf(&text[intoffset], "%14lg", output) != 1) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                  "Was expecting to read a double from file line: %s", text);
      return FALSE;
    }
  }
  return TRUE;
}

/** Read in a single MDT feature definition from the file. Return TRUE
    on success. */
static gboolean read_mdt_feature(struct mdt *mdt,
                                 const struct mdt_library *mlib, int nfeat,
                                 int ifeat, int nbins, GError **err)
{
  struct mod_mdt_feature *feat;

  mod_mdt_nfeat_set(&mdt->base, nfeat);
  feat = &mdt->base.features[nfeat - 1];
  feat->ifeat = ifeat;
  feat->istart = 1;
  feat->nbins = feat->iend = nbins;

  if (!check_feature_type(feat->ifeat, mlib, err)) {
    return FALSE;
  } else {
    struct mod_mdt_libfeature *libfeat = &mlib->base.features[ifeat - 1];
    if (nbins != libfeat->nbins) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                  "MDT number of bins (%d) is different to that in the bin "
                  "file (%d) for feature number %d, feature type %d", nbins,
                  libfeat->nbins, nfeat - 1, feat->ifeat);
      return FALSE;
    }
  }
  return TRUE;
}

/** Read the definition of MDT features from the file. Return TRUE on
    success. */
static gboolean read_mdt_features(struct mdt *mdt,
                                  const struct mdt_library *mlib, FILE *fp,
                                  GError **err)
{
  int nfeat = 0, eof;
  gboolean retval = TRUE;
  GString *str = g_string_new(NULL);

  mod_mdt_nfeat_set(&mdt->base, nfeat);
  while (retval) {
    int ifeat, nbins;

    if (!mdt_file_read_line(fp, str, &eof, err)) {
      retval = FALSE;
    } else if (eof != 0 || strlen(str->str) == 0) {
      break;
    } else if (strlen(str->str) > 6
               && sscanf(&str->str[6], "%6d%6d", &ifeat, &nbins) == 2) {
      retval = read_mdt_feature(mdt, mlib, ++nfeat, ifeat, nbins, err);
    } else {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                  "Was expecting to read feature information from "
                  "MDT line: %s", str->str);
      retval = FALSE;
    }
  }
  g_string_free(str, TRUE);
  return retval;
}

/** Read the MDT feature offsets and shape from the file. Return TRUE on
    success. */
static gboolean read_mdt_offsets(struct mdt *mdt, FILE *fp, GError **err)
{
  int i;
  gboolean retval = TRUE;
  GString *str = g_string_new(NULL);

  for (i = 0; i < mdt->base.nfeat && retval; i++) {
    int eof, nfeat, istart, iend;
    if (!mdt_file_read_line(fp, str, &eof, err)) {
      retval = FALSE;
    } else if (eof != 0) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                  "End of file encountered while reading feature offsets "
                  "and shape: expected to read data for %d features",
                  mdt->base.nfeat);
      retval = FALSE;
    } else if (sscanf(str->str, "%3d %6d %6d", &nfeat, &istart, &iend) != 3) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                  "Was expecting to read feature offset and shape: got %s",
                  str->str);
      retval = FALSE;
    } else {
      struct mod_mdt_feature *feat = &mdt->base.features[i];
      feat->istart = istart;
      feat->iend = iend;
    }
  }

  g_string_free(str, TRUE);
  return retval;
}

/** Read a single line of MDT data from the file. Return TRUE on success. */
static gboolean read_mdt_data_line(FILE *fp, GString *str, double *output,
                                   GError **err)
{
  int eof;

  if (!mdt_file_read_line(fp, str, &eof, err)) {
    return FALSE;
  } else if (eof != 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                "End of file while reading MDT bin data");
    return FALSE;
  } else if (sscanf(str->str, "%15lg", output) != 1) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                "Expected to read a double: got %s", str->str);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Read the footer from the MDT file. Return TRUE on success. */
static gboolean read_mdt_footer(FILE *fp, GString *str, gboolean *pdf,
                                GError **err)
{
  int eof;

  if (!mdt_file_read_line(fp, str, &eof, err)) {
    return FALSE;
  } else if (eof != 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                "End of file while reading MDT footer");
    return FALSE;
  } else if (strcmp(str->str, "MDT TABLE END:PDF") == 0) {
    *pdf = TRUE;
    return TRUE;
  } else if (strcmp(str->str, "MDT TABLE END") == 0) {
    *pdf = FALSE;
    return TRUE;
  } else {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT, "Corrupt MDT footer");
    return FALSE;
  }
}

/** Read the MDT data from the file. Return TRUE on success. */
static gboolean read_mdt_data(struct mdt *mdt, FILE *fp, const char *header,
                              GError **err)
{
  int nelems;
  if (strlen(header) > 16 && sscanf(&header[16], "%9d", &nelems) == 1) {
    int i;
    gboolean retval = TRUE;
    GString *str = g_string_new(NULL);

    mod_mdt_nelems_set(&mdt->base, nelems);
    for (i = 0; i < nelems && retval; i++) {
      retval = read_mdt_data_line(fp, str, &mdt->base.bin[i], err);
    }

    if (retval) {
      retval = read_mdt_footer(fp, str, &mdt->pdf, err);
    }
    g_string_free(str, TRUE);
    return retval;
  } else {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FILE_FORMAT,
                "Expecting to read MDT table number of elements: got %s",
                header);
    return FALSE;
  }
}

/** Read a single line of text from an MDT file. Return TRUE on success. */
static gboolean read_mdt_file_line(struct mdt *mdt,
                                   const struct mdt_library *mlib, FILE *fp,
                                   GString *str, GError **err)
{
  if (strcmp(str->str, "  #  FEATURE NBINS NAME") == 0) {
    return read_mdt_features(mdt, mlib, fp, err);
  } else if (strcmp(str->str, "  # ISTART   IEND") == 0) {
    return read_mdt_offsets(mdt, fp, err);
  } else if (strncmp(str->str, "MDT TABLE START:", 16) == 0) {
    return read_mdt_data(mdt, fp, str->str, err);
  } else if (!read_int_line(str->str, "Number of alignments", 39,
                            &mdt->nalns, err)
             || !read_int_line(str->str,
                               "Number of proteins in the alignments", 39,
                               &mdt->n_proteins, err)
             || !read_int_line(str->str, "Number of proteins or their pairs",
                               39, &mdt->n_protein_pairs, err)
             || !read_double_line(str->str, "Sample size", 39,
                                  &mdt->sample_size, err)) {
    return FALSE;
  }
  return TRUE;
}

/** Read in an MDT in text format. Return TRUE on success. */
gboolean mdt_read(struct mdt *mdt, const struct mdt_library *mlib,
                  const char *filename, GError **err)
{
  struct mod_file file_info;
  FILE *fp;

  /* First, blank out any existing features */
  mod_mdt_nfeat_set(&mdt->base, 0);

  fp = mdt_open_file(filename, "r", &file_info, err);
  if (fp) {
    int eof = 0;
    gboolean retval = TRUE;
    GString *str = g_string_new(NULL);

    while (eof == 0 && retval) {
      if (!mdt_file_read_line(fp, str, &eof, err)) {
        retval = FALSE;
      } else if (eof == 0) {
        retval = read_mdt_file_line(mdt, mlib, fp, str, err);
      }
    }

    g_string_free(str, TRUE);
    retval &= mdt_close_file(fp, &file_info, err);

    if (retval) {
      retval = mdt_setup(mdt, mlib, err);
    }
    return retval;
  } else {
    return FALSE;
  }
}
