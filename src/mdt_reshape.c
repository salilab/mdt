/** \file mdt_reshape.c    Functions to reshape MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Set new MDT indices using new offset and shape */
static void reshape_mdt_indices(const struct mdt_type *mdtin,
                                struct mdt_type *mdtout, const int offset[],
                                const int shape[], const int features[],
                                const int old_position[])
{
  int i;
  int *in_iend = f_int1_pt(&mdtin->iend);
  int *istart = f_int1_pt(&mdtout->istart);
  int *iend = f_int1_pt(&mdtout->iend);
  int *ifeat = f_int1_pt(&mdtout->ifeat);
  int *nbins = f_int1_pt(&mdtout->nbins);

  mdtout->nelems = 1;
  for (i = 0; i < mdtin->nfeat; i++) {
    istart[i] = offset[i] + 1;
    if (shape[i] <= 0) {
      iend[i] = in_iend[old_position[i]] + shape[i];
    } else {
      iend[i] = offset[i] + shape[i];
    }

    ifeat[i] = features[i];
    nbins[i] = iend[i] - istart[i] + 1;
    mdtout->nelems *= nbins[i];
  }

  make_mdt_stride(mdtout);
}


/** Do the hard work of reshaping the table. */
static void reshape_mdt_table(const struct mdt_type *mdtin,
                              struct mdt_type *mdtout, const int new_position[])
{
  int *out_indf, *in_indf;
  double *out_bin, *in_bin;
  out_bin = f_double1_pt(&mdtout->bin);
  in_bin = f_double1_pt(&mdtin->bin);
  out_indf = mdt_start_indices(mdtout);
  in_indf = malloc(sizeof(int) * mdtin->nfeat);
  do {
    int i, i1, i2;
    for (i = 0; i < mdtin->nfeat; i++) {
      in_indf[i] = out_indf[new_position[i]];
    }
    i1 = indmdt(in_indf, mdtin);
    i2 = indmdt(out_indf, mdtout);
    out_bin[i2] = in_bin[i1];
  } while (roll_ind(out_indf, f_int1_pt(&mdtout->istart),
                    f_int1_pt(&mdtout->iend), mdtout->nfeat));
  free(in_indf);
  free(out_indf);
}

/** Get mapping from old to new features */
static void get_position_mappings(const struct mdt_type *mdt,
                                  const int features[], int old_position[],
                                  int new_position[], const char *routine,
                                  int *ierr)
{
  int i, j, *ifeat;
  ifeat = f_int1_pt(&mdt->ifeat);
  for (i = 0; i < mdt->nfeat; i++) {
    int match = 0;
    for (j = 0; j < mdt->nfeat; j++) {
      if (features[i] == ifeat[j]) {
        old_position[i] = j;
        new_position[j] = i;
        match = 1;
        break;
      }
    }
    if (!match) {
      modlogerror(routine, ME_VALUE,
                  "Feature type %d does not exist in input MDT.", features[i]);
      *ierr = 1;
      return;
    }
  }
}
 

/** Check new offset and shape */
static void check_start_end(const struct mdt_type *mdt, const int offset[],
                            const int shape[], const int old_position[],
                            const int features[], const char *routine,
                            int *ierr)
{
  int i, *iend, *istart, *nbins;
  iend = f_int1_pt(&mdt->iend);
  istart = f_int1_pt(&mdt->istart);
  nbins = f_int1_pt(&mdt->nbins);
  for (i = 0; i < mdt->nfeat; i++) {
    int end;
    if (shape[i] <= 0) {
      end = iend[old_position[i]] + shape[i];
    } else {
      end = offset[i] + shape[i];
    }
    if (offset[i] + 1 < istart[old_position[i]] || end < offset[i] + 1
        || end > nbins[old_position[i]]) {
      modlogerror(routine, ME_GENERIC,
                  "For feature %d, new start %d and size %d are out of range.",
                  features[i], offset[i], shape[i]);
      *ierr = 1;
      return;
    }
  }
}


/** Reshape an MDT. */
void mdt_reshape(const struct mdt_type *mdtin, struct mdt_type *mdtout,
                 const int features[], int n_features, const int offset[],
                 int n_offset, const int shape[], int n_shape, int *ierr)
{
  const char *routine = "mdt_reshape";
  int *old_position, *new_position;

  *ierr = 0;
  if (n_features != mdtin->nfeat || n_offset != mdtin->nfeat
      || n_shape != mdtin->nfeat) {
    modlogerror(routine, ME_VALUE, "features, offset and shape must all match"
                " the dimension of the MDT (%d)", mdtin->nfeat);
    *ierr = 1;
    return;
  }

  old_position = malloc(sizeof(int) * mdtin->nfeat);
  new_position = malloc(sizeof(int) * mdtin->nfeat);

  get_position_mappings(mdtin, features, old_position, new_position, routine,
                        ierr);
  if (*ierr) {
    free(old_position);
    free(new_position);
    return;
  }
  check_start_end(mdtin, offset, shape, old_position, features, routine, ierr);
  if (*ierr) {
    free(old_position);
    free(new_position);
    return;
  }

  copy_mdt(mdtin, mdtout);
  reshape_mdt_indices(mdtin, mdtout, offset, shape, features, old_position);

  /* reshape the MDT table: */
  reshape_mdt_table(mdtin, mdtout, new_position);

  /* a little heuristic here: */
  if (!mdtout->pdf) {
    mdtout->sample_size = get_sum(f_double1_pt(&mdtout->bin), mdtout->nelems);
  }

  free(old_position);
  free(new_position);
}
