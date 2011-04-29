/** \file atom_hydrogen_bonds.c  Atom hydrogen bond features.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../geometry.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"
#include "../mdt_hydrogen_bonds.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  int hbprop_type = GPOINTER_TO_INT(feat->data);
  const int *binprop;
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float1_get(&s->cd.x, atom);
  if (coordinate_undefined(f)) {
    return mdt_feature_undefined_bin_get(feat);
  } else if ((binprop = property_hb_iatta(aln, protein, prop, mlib,
                                          libs, err))) {
    return feat_to_bin(numb_hda(atom, binprop, &s->cd, mlib->hbond,
                                mlib->hbond_cutoff, hbprop_type), feat);
  } else {
    return -1;
  }
}

static int make_feature(struct mdt_library *mlib, gboolean pos2,
                        const char *name, int hbprop_type)
{
  return mdt_feature_atom_add(mlib, name, MOD_MDTC_NONE, pos2, getbin,
                              GINT_TO_POINTER(hbprop_type), NULL);
}

int mdt_feature_hydrogen_bond_donor(struct mdt_library *mlib, gboolean pos2)
{
  return make_feature(mlib, pos2, "H-bond donor", 0);
}

int mdt_feature_hydrogen_bond_acceptor(struct mdt_library *mlib, gboolean pos2)
{
  return make_feature(mlib, pos2, "H-bond acceptor", 1);
}

int mdt_feature_hydrogen_bond_charge(struct mdt_library *mlib, gboolean pos2)
{
  return make_feature(mlib, pos2, "Total H-bond charge around atom", 2);
}
