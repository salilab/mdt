/** \file tuple_angle.c     Tuple-tuple non-bonded angle features.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"
#include "../mdt_atom_classes.h"
#include "../mdt_tuples.h"
#include "../geometry.h"

static int angle1(const struct mod_alignment *aln, int protein,
                  int atom1, const struct mdt_tuple *tuple1,
                  int atom2, const struct mdt_tuple *tuple2,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return iangle0(atom1, atom2, tuple2->iata[0], s, feat);
}

static int angle2(const struct mod_alignment *aln, int protein,
                  int atom1, const struct mdt_tuple *tuple1,
                  int atom2, const struct mdt_tuple *tuple2,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return iangle0(tuple1->iata[0], atom1, atom2, s, feat);
}

static int angle_feature(struct mdt_library *mlib, const char *name,
                         mdt_cb_feature_tuple_pair getbin, GError **err)
{
  int ifeat;
  if (!tuple_require_natom(mlib, 2,  err)) {
    return -1;
  }
  ifeat = mdt_feature_tuple_pair_add(mlib, name, MOD_MDTC_NONE, getbin, NULL,
                                     NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}

int mdt_feature_tuple_angle1(struct mdt_library *mlib, GError **err)
{
  return angle_feature(mlib, "Tuple-tuple non-bonded angle1", angle1, err);
}

int mdt_feature_tuple_angle2(struct mdt_library *mlib, GError **err)
{
  return angle_feature(mlib, "Tuple-tuple non-bonded angle2", angle2, err);
}
