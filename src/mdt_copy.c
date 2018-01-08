/** \file mdt_copy.c       Functions to copy MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2018 Andrej Sali
 */

#include <stdlib.h>
#include <string.h>
#include "modeller.h"
#include "mdt.h"

/** Copy the base Modeller MDT */
static void mod_mdt_copy(const struct mod_mdt *mdtin, struct mod_mdt *mdtout,
                         mod_mdt_bin_type bin_type)
{
  int i;
  mdtout->nbinx = mdtin->nbinx;
  mdtout->nfeat1 = mdtin->nfeat1;
  mdtout->bin_type = bin_type;
  for (i = 0; i < MAXDATATYP; i++) {
    mdtout->readin[i] = mdtin->readin[i];
  }
  mdtout->nprotcmp = mdtin->nprotcmp;

  mod_mdt_nfeat_set(mdtout, mdtin->nfeat);
  mod_mdt_nelems_set(mdtout, mdtin->nelems);
  for (i = 0; i < mdtin->nfeat; i++) {
    struct mod_mdt_feature *in, *out;
    in = &mdtin->features[i];
    out = &mdtout->features[i];
    out->stride = in->stride;
    out->ifeat = in->ifeat;
    out->istart = in->istart;
    out->iend = in->iend;
    out->nbins = in->nbins;
  }
  /* If types are the same, copy bin data straight across */
  if (mdtout->bin_type == mdtin->bin_type) {
    memcpy(mdtout->bindata, mdtin->bindata,
           mod_mdt_bin_get_size(mdtin) * mdtin->nelems);
  } else {
    /* Otherwise, convert each data point to the new storage type */
    for (i = 0; i < mdtin->nelems; ++i) {
      mod_mdt_bin_set(mdtout, i, mod_mdt_bin_get(mdtin, i));
    }
  }
}

static void insert_func(gpointer key, gpointer value, gpointer user_data)
{
  GHashTable *t = (GHashTable *)user_data;
  g_hash_table_insert(t, key, value);
}

/** Make mdtout a copy of mdtin. */
void mdt_copy(const struct mdt *mdtin, struct mdt *mdtout,
              mod_mdt_bin_type bin_type)
{
  mod_mdt_copy(&mdtin->base, &mdtout->base, bin_type);
  mdtout->scantype = mdtin->scantype;
  mdtout->symmetric = mdtin->symmetric;
  mdtout->nalns = mdtin->nalns;
  mdtout->n_protein_pairs = mdtin->n_protein_pairs;
  mdtout->n_proteins = mdtin->n_proteins;
  mdtout->sample_size = mdtin->sample_size;
  mdtout->pdf = mdtin->pdf;
  memcpy(&mdtout->scan_params, &mdtin->scan_params,
         sizeof(struct mdt_scan_parameters));
  g_hash_table_foreach(mdtin->write_lib_funcs, insert_func,
                       mdtout->write_lib_funcs);
}
