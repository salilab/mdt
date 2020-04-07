/** \file bond_type.c  Bond type feature.
 *
 *             Part of MDT, Copyright(c) 1989-2020 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_atom_classes.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  const struct mdt_bond *bond, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  return bond->bndgrp;
}

static int add_feature(struct mdt_library *mlib, int ifeat)
{
  struct mdt_feature *feat = &g_array_index(mlib->features, struct mdt_feature,
                                            ifeat - 1);
  struct mod_mdt_libfeature *libfeat = &mlib->base.features[ifeat - 1];

  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  /* Set number of bins and their symbols */
  update_mdt_feat_atclass(libfeat, mlib->atclass[feat->u.bond.type + 1]);
  return ifeat;
}

int mdt_feature_bond_type(struct mdt_library *mlib)
{
  return add_feature(mlib, mdt_feature_bond_add(mlib, "Bond type",
                                                MOD_MDTC_NONE, getbin, NULL,
                                                NULL));
}

int mdt_feature_angle_type(struct mdt_library *mlib)
{
  return add_feature(mlib, mdt_feature_angle_add(mlib, "Angle type",
                                                 MOD_MDTC_NONE, getbin, NULL,
                                                 NULL));
}

int mdt_feature_dihedral_type(struct mdt_library *mlib)
{
  return add_feature(mlib,
                     mdt_feature_dihedral_add(mlib, "Dihedral type",
                                              MOD_MDTC_NONE, getbin, NULL,
                                              NULL));
}
