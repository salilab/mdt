/** \file gap_distance.c  Distance to a gap features.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

/** Get distance to a gap from the given alignment position,
    in the given direction. */
static int gap_distance_dir(const struct mod_alignment *aln, int protein1,
                            int protein2, int alnpos, int inc)
{
  int i, gapdist = 0;
  for (i = alnpos; i >= 0 && i < aln->naln; i += inc) {
    int ir1 = mod_int2_get(&aln->ialn, i, protein1);
    int ir2 = mod_int2_get(&aln->ialn, i, protein2);
    /* If residues at *both* positions, increase distance to gap */
    if (ir1 > 0 && ir2 > 0) {
      gapdist++;
    /* Otherwise, if residue at only one position, we've found the gap */
    } else if (ir1 > 0 || ir2 > 0) {
      break;
    }
    /* Otherwise, there are gaps in both sequences. Skip over it. */
  }
  return gapdist;
}

/** Get distance to a gap from the given alignment position. */
static int gap_distance(const struct mod_alignment *aln, int protein1,
                        int protein2, int alnpos)
{
  return MIN(gap_distance_dir(aln, protein1, protein2, alnpos, 1),
             gap_distance_dir(aln, protein1, protein2, alnpos, -1));
}

static int getbin(const struct mod_alignment *aln, int protein1, int protein2,
                  int alnpos, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float f = gap_distance(aln, protein1, protein2, alnpos);
  return feat_to_bin(f, feat);
}

static int avgetbin(const struct mod_alignment *aln, int protein1, int protein2,
                    int alnpos1, int alnpos2, struct mdt_properties *prop,
                    const struct mdt_feature *feat,
                    const struct mdt_library *mlib,
                    const struct mod_libraries *libs, GError **err)
{
  float f1 = gap_distance(aln, protein1, protein2, alnpos1);
  float f2 = gap_distance(aln, protein1, protein2, alnpos2);
  return feat_to_bin(0.5 * (f1 + f2), feat);
}

int mdt_feature_gap_distance(struct mdt_library *mlib, int protein1,
                             int protein2, GError **err)
{
  return mdt_feature_aligned_residue_add(mlib, "Distance from a gap",
                                         MOD_MDTC_NONE, protein1, protein2,
                                         getbin, NULL, NULL, err);
}

int mdt_feature_average_gap_distance(struct mdt_library *mlib, int protein1,
                                     int protein2, GError **err)
{
  return mdt_feature_aligned_residue_pair_add(mlib,
                                              "Average distance from a gap",
                                              MOD_MDTC_NONE, protein1, protein2,
                                              FALSE, avgetbin, NULL, NULL, err);
}
