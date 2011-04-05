/** \file tuple_dihedral.c     Tuple-tuple non-bonded dihedral features.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"
#include "../mdt_atom_classes.h"
#include "../mdt_tuples.h"
#include "../geometry.h"

static int dihedral1(const struct mod_alignment *aln, int protein,
                     int atom1, const struct mdt_tuple *tuple1,
                     int atom2, const struct mdt_tuple *tuple2,
                     struct mdt_properties *prop,
                     const struct mdt_feature *feat,
                     const struct mdt_library *mlib,
                     const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idihedral0(tuple1->iata[0], atom1, atom2, tuple2->iata[0], s, feat);
}

static int dihedral2(const struct mod_alignment *aln, int protein,
                     int atom1, const struct mdt_tuple *tuple1,
                     int atom2, const struct mdt_tuple *tuple2,
                     struct mdt_properties *prop,
                     const struct mdt_feature *feat,
                     const struct mdt_library *mlib,
                     const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idihedral0(tuple1->iata[1], tuple2->iata[0], atom1, atom2, s, feat);
}

static int dihedral3(const struct mod_alignment *aln, int protein,
                     int atom1, const struct mdt_tuple *tuple1,
                     int atom2, const struct mdt_tuple *tuple2,
                     struct mdt_properties *prop,
                     const struct mdt_feature *feat,
                     const struct mdt_library *mlib,
                     const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idihedral0(atom1, atom2, tuple1->iata[0], tuple2->iata[1], s, feat);
}

static int dihedral_feature(struct mdt_library *mlib, const char *name,
                            mdt_cb_feature_tuple_pair getbin, int need_natom,
                            GError **err)
{
  int ifeat;
  if (!tuple_require_natom(mlib, need_natom, err)) {
    return -1;
  }
  ifeat = mdt_feature_tuple_pair_add(mlib, name, MOD_MDTC_NONE, getbin, NULL,
                                     NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_periodic_set(mlib, ifeat, TRUE);
  return ifeat;
}

int mdt_feature_tuple_dihedral1(struct mdt_library *mlib, GError **err)
{
  return dihedral_feature(mlib, "Tuple-tuple non-bonded dihedral1", dihedral1,
                          2, err);
}

int mdt_feature_tuple_dihedral2(struct mdt_library *mlib, GError **err)
{
  return dihedral_feature(mlib, "Tuple-tuple non-bonded dihedral2", dihedral2,
                          3, err);
}

int mdt_feature_tuple_dihedral3(struct mdt_library *mlib, GError **err)
{
  return dihedral_feature(mlib, "Tuple-tuple non-bonded dihedral3", dihedral3,
                          3, err);
}
