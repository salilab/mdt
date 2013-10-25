/** \file mdt_types.c  Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_index.h"
#include "mdt_types.h"
#include "mdt_atom_classes.h"
#include "mdt_feature.h"
#include "mdt_residue_bonds.h"

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
  mdt->write_lib_funcs = g_hash_table_new(NULL, NULL);
  mdt->scan_params.scan_called = FALSE;
  return mdt;
}

void mdt_set_write_lib_callback(struct mdt *mdt, mdt_cb_write_lib writelibfunc)
{
  g_hash_table_insert(mdt->write_lib_funcs, writelibfunc, GINT_TO_POINTER(1));
}

/** Free an mdt structure */
void mdt_free(struct mdt *mdt)
{
  mod_mdt_dealloc(&mdt->base);
  g_hash_table_destroy(mdt->write_lib_funcs);
  g_free(mdt);
}

/** Make a new mdt_library structure */
struct mdt_library *mdt_library_new(struct mod_libraries *libs)
{
  int i;
  struct mdt_library *mlib;
  mlib = g_malloc(sizeof(struct mdt_library));
  mod_mdt_library_init(&mlib->base);
  mlib->libs = libs;
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
  mlib->distance_atoms[0] = g_strdup("");
  mlib->distance_atoms[1] = g_strdup("");
  mdt_residue_bond_list_init(&mlib->residue_bond_list);
  mlib->atom_properties = g_array_new(FALSE, FALSE,
                                      sizeof(struct mdt_user_property));
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
  for (i = 0; i < mlib->base.nfeat; i++) {
    struct mdt_feature *feat;
    feat = &g_array_index(mlib->features, struct mdt_feature, i);
    if (feat->freefunc) {
      feat->freefunc(feat->data);
    }
  }
  mdt_atom_class_list_free(mlib->hbond);
  mdt_atom_class_list_free(mlib->tupclass);
  g_array_free(mlib->features, TRUE);
  g_free(mlib->distance_atoms[0]);
  g_free(mlib->distance_atoms[1]);
  mdt_residue_bond_list_free(&mlib->residue_bond_list);
  g_array_free(mlib->atom_properties, TRUE);
  g_free(mlib);
}

int mdt_library_add_atom_property(struct mdt_library *mlib,
                                  mdt_cb_get_property get_property,
                                  gpointer data, GDestroyNotify freefunc)
{
  struct mdt_user_property p;
  p.get_property = get_property;
  p.data = data;
  p.freefunc = freefunc;
  g_array_append_val(mlib->atom_properties, p);
  return mlib->atom_properties->len - 1;
}
