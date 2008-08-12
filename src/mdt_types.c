/** \file mdt_types.c  Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_index.h"
#include "mdt_types.h"
#include "mdt_atom_classes.h"
#include "mdt_feature.h"

/** Make a new mdt structure */
struct mdt *mdt_new(mod_mdt_bin_type bin_type)
{
  struct mdt *mdt;
  mdt = g_malloc(sizeof(struct mdt));
  mod_mdt_init(&mdt->base);
  mdt->pdf = FALSE;
  mdt->sample_size = 0.;
  mdt->nalns = mdt->n_proteins = mdt->n_protein_pairs = 0;
  mdt->symmetric = FALSE;
  mdt->scantype = 0;
  mdt->base.bin_type = bin_type;
  return mdt;
}

/** Free an mdt structure */
void mdt_free(struct mdt *mdt)
{
  mod_mdt_dealloc(&mdt->base);
  g_free(mdt);
}

/** Make a new mdt_library structure */
struct mdt_library *mdt_library_new(void)
{
  int i;
  struct mdt_library *mlib;
  mlib = g_malloc(sizeof(struct mdt_library));
  mod_mdt_library_init(&mlib->base);
  mlib->hbond_cutoff = 3.5;
  mlib->special_atoms = FALSE;
  mlib->feature_added = FALSE;
  for (i = 0; i < 4; i++) {
    mlib->atclass[i] = mdt_atom_class_list_new(i + 1);
  }
  mlib->hbond = mdt_atom_class_list_new(1);

  /* Set natom to 0 to start with; gets set to 2 or 3 when we read the
     doublet/triplet class file */
  mlib->tupclass = mdt_atom_class_list_new(0);
  mlib->features = g_array_new(FALSE, TRUE, sizeof(struct mdt_feature));
  mlib->features = g_array_set_size(mlib->features, mlib->base.nfeat);
  return mlib;
}

/** Free an mdt_library structure */
void mdt_library_free(struct mdt_library *mlib)
{
  int i;
  mod_mdt_library_dealloc(&mlib->base);
  for (i = 0; i < 4; i++) {
    mdt_atom_class_list_free(mlib->atclass[i]);
  }
  mdt_atom_class_list_free(mlib->hbond);
  mdt_atom_class_list_free(mlib->tupclass);
  g_array_free(mlib->features, TRUE);
  g_free(mlib);
}
