/** \file mdt_read_hdf5.c   Functions to read MDTs from binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_hdf5.h"
#include "mdt.h"
#include "util.h"

/** Read MDT feature information from an HDF5 file. Return TRUE on success. */
static gboolean read_mdt_features(hid_t file_id, struct mod_mdt *mdt,
                                  GError **err)
{
  hsize_t nfeat, *shape;
  int *ifeat, *offset;
  herr_t ret;

  ret = mod_dataset_get_ndsize(file_id, "/features", 1, &nfeat);
  if (ret < 0) {
    handle_modeller_hdf5_error(err);
    return FALSE;
  }
  mod_mdt_nfeat_set(mdt, nfeat);
  ifeat = g_malloc(mdt->nfeat * sizeof(int));
  offset = g_malloc(mdt->nfeat * sizeof(int));
  shape = g_malloc(mdt->nfeat * sizeof(hsize_t));
  if (mod_dataset_read_int(file_id, "/features", 1, &nfeat, ifeat) >= 0
      && mod_dataset_read_int(file_id, "/offset", 1, &nfeat, offset) >= 0
      && mod_dataset_get_ndsize(file_id, "/mdt", nfeat, shape) >= 0) {
    int i;
    for (i = 0; i < mdt->nfeat; i++) {
      struct mod_mdt_feature *feat = &mdt->features[i];
      feat->ifeat = ifeat[i];
      feat->istart = offset[i] + 1;
      feat->nbins = shape[i];
      feat->iend = feat->istart + feat->nbins - 1;
    }
  } else {
    handle_modeller_hdf5_error(err);
    ret = -1;
  }
  g_free(ifeat);
  g_free(offset);
  g_free(shape);
  return (ret >= 0);
}


/** Check a single MDT feature for sanity. */
static gboolean check_mdt_feature(const struct mod_mdt_feature *feat,
                                  int nbins, int i,
                                  const struct mdt_library *mlib, GError **err)
{
  const struct mod_mdt_libfeature *libfeat;

  if (!check_feature_type(feat->ifeat, mlib, err)) {
    return FALSE;
  }

  libfeat = &mlib->base.features[feat->ifeat - 1];
  if (nbins != libfeat->nbins) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "MDT number of bins (%d) is different to that in the bin "
                "file (%d) for feature number %d, feature type %d", nbins,
                libfeat->nbins, i, feat->ifeat);
    return FALSE;
  } else if (feat->istart < 1 || feat->istart > libfeat->nbins) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "MDT table offset (%d) is outside of bin file range 0-%d "
                "for feature number %d, feature type %d", feat->istart - 1,
                libfeat->nbins - 1, i, feat->ifeat);
    return FALSE;
  } else if (feat->iend < 1 || feat->iend > libfeat->nbins) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "MDT table end (%d) is outside of bin file range 0-%d "
                "for feature number %d, feature type %d", feat->iend - 1,
                libfeat->nbins - 1, i, feat->ifeat);
    return FALSE;
  } else {
    return TRUE;
  }
}


/** Check all MDT feature information for sanity. */
static gboolean check_mdt_features(hid_t file_id, const struct mod_mdt *mdt,
                                   const struct mdt_library *mlib,
                                   GError **err)
{
  int i, *nbins;
  gboolean retval;
  hsize_t nfeat = mdt->nfeat;

  nbins = g_malloc(mdt->nfeat * sizeof(int));

  retval = (mod_dataset_read_int(file_id, "/nbins", 1, &nfeat, nbins) >= 0);
  if (!retval) {
    handle_modeller_hdf5_error(err);
  }
  for (i = 0; i < mdt->nfeat && retval; i++) {
    retval = check_mdt_feature(&mdt->features[i], nbins[i], i, mlib, err);
  }

  g_free(nbins);
  return retval;
}


/** Read MDT data from an HDF5 file. Return TRUE on success. */
static gboolean read_mdt_data(hid_t file_id, struct mdt *mdt, GError **err)
{
  herr_t ret;
  hsize_t *dims;
  int i, nelems;

  dims = g_malloc(mdt->base.nfeat * sizeof(hsize_t));
  nelems = 1;
  for (i = 0; i < mdt->base.nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->base.features[i];
    dims[i] = feat->nbins;
    nelems *= feat->nbins;
  }

  mod_mdt_nelems_set(&mdt->base, nelems);
  if (mdt->base.bin_type == MOD_MDTB_FLOAT) {
    ret = mod_dataset_read_float(file_id, "/mdt", mdt->base.nfeat, dims,
                                 (float *)mdt->base.bindata);
  } else {
    ret = mod_dataset_read_double(file_id, "/mdt", mdt->base.nfeat, dims,
                                  (double *)mdt->base.bindata);
  }
  if (ret >= 0) {
    hsize_t featdim = 1;
    char is_pdf;
    if (mod_dataset_read_int(file_id, "/n_alignments", 1, &featdim,
                             &mdt->nalns) < 0
        || mod_dataset_read_int(file_id, "/n_proteins", 1, &featdim,
                                &mdt->n_proteins) < 0
        || mod_dataset_read_int(file_id, "/n_protein_pairs", 1, &featdim,
                                &mdt->n_protein_pairs) < 0
        || mod_dataset_read_double(file_id, "/sample_size", 1, &featdim,
                                   &mdt->sample_size) < 0
        || mod_dataset_read_char(file_id, "/pdf", 1, &featdim, &is_pdf) < 0) {
      ret = -1;
    } else {
      mdt->pdf = (is_pdf == 1);
    }
  }
  g_free(dims);

  if (ret < 0) {
    handle_modeller_hdf5_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Read in an MDT in HDF5 format. Return TRUE on success. */
gboolean mdt_read_hdf5(struct mdt *mdt, const struct mdt_library *mlib,
                       const char *filename, GError **err)
{
  hid_t file_id;
  struct mod_file file_info;

  file_id = mdt_hdf_open(filename, H5F_ACC_RDONLY, H5P_DEFAULT, &file_info,
                         err);
  if (file_id < 0 || !read_mdt_features(file_id, &mdt->base, err)
      || !check_mdt_features(file_id, &mdt->base, mlib, err)
      || !read_mdt_data(file_id, mdt, err)
      || !mdt_hdf_close(file_id, &file_info, err)) {
    return FALSE;
  } else {
    return mdt_setup(mdt, mlib, err);
  }
}
