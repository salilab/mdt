/** \file mdt_copy.c       Functions to copy MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include <string.h>
#include "modeller.h"
#include "mdt.h"

/** Copy the base Modeller MDT */
static void mod_mdt_copy(const struct mod_mdt *mdtin, struct mod_mdt *mdtout)
{
  int i;
  mdtout->nalns = mdtin->nalns;
  mdtout->n_protein_pairs = mdtin->n_protein_pairs;
  mdtout->n_proteins = mdtin->n_proteins;
  mdtout->nbinx = mdtin->nbinx;
  mdtout->sample_size = mdtin->sample_size;
  mdtout->nfeat1 = mdtin->nfeat1;
  mdtout->pdf = mdtin->pdf;
  for (i = 0; i < MAXDATATYP; i++) {
    mdtout->readin[i] = mdtin->readin[i];
  }
  mdtout->symmetric = mdtin->symmetric;
  mdtout->nprotcmp = mdtin->nprotcmp;
  mdtout->nresfeat = mdtin->nresfeat;

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
  memcpy(mdtout->bin, mdtin->bin, sizeof(double) * mdtin->nelems);
}

/** Make mdtout a copy of mdtin. */
void mdt_copy(const struct mdt *mdtin, struct mdt *mdtout)
{
  mod_mdt_copy(&mdtin->base, &mdtout->base);
}
