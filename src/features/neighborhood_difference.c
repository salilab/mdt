/** \file neighborhood_difference.c  Neighborhood difference features.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

/** Get the neighborhood difference for the given alignment position. This is
    simply the average of the residue distance scores for all neighbors. */
static float local_difference(const struct mod_alignment *aln, int protein1,
                              int protein2, int alnpos,
                              const struct mod_libraries *libs)
{
  struct mod_structure *s1 = mod_alignment_structure_get(aln, protein1);
  struct mod_sequence *seq1 = mod_alignment_sequence_get(aln, protein1);
  struct mod_sequence *seq2 = mod_alignment_sequence_get(aln, protein2);
  int ires = mod_int2_get(&aln->ialn, alnpos, protein1) - 1;
  int neigh = mod_int1_get(&s1->neigh, ires);
  float score = 0.;
  int i;
  for (i = 0; i < neigh; ++i) {
    /* Map each neighbor residue onto the other sequence */
    int ires1 = mod_int2_get(&s1->ineigh, i, ires) - 1;
    int ipos = mod_int2_get(&aln->invaln, ires1, protein1) - 1;
    int ires2 = mod_int2_get(&aln->ialn, ipos, protein2) - 1;
    /* Get types of the residues in both sequences, and map
       any nonstandard residues to types valid for the RR matrix */
    int irtyp1 = mod_int1_get(&seq1->irestyp, ires1);
    int irtyp2;
    if (ires2 < 0) {
      irtyp2 = libs->igaptyp;
    } else {
      irtyp2 = mod_int1_get(&seq2->irestyp, ires2);
      irtyp2 = mod_residue_standard_type(irtyp2, libs);
    }
    irtyp1 = mod_residue_standard_type(irtyp1, libs);
    /* Get difference score for this residue pair */
    score += mod_float2_get(&libs->rrwght, irtyp1 - 1, irtyp2 - 1);
  }
  if (neigh > 0) {
    score /= neigh;
  }
  return score;
}

static int getbin(const struct mod_alignment *aln, int protein1, int protein2,
                  int alnpos, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float f = local_difference(aln, protein1, protein2, alnpos, libs);
  return feat_to_bin(f, feat);
}

static int avgetbin(const struct mod_alignment *aln, int protein1, int protein2,
                    int alnpos1, int alnpos2, struct mdt_properties *prop,
                    const struct mdt_feature *feat,
                    const struct mdt_library *mlib,
                    const struct mod_libraries *libs, GError **err)
{
  float f1 = local_difference(aln, protein1, protein2, alnpos1, libs);
  float f2 = local_difference(aln, protein1, protein2, alnpos2, libs);
  return feat_to_bin(0.5 * (f1 + f2), feat);
}

int mdt_feature_neighborhood_difference(struct mdt_library *mlib, int protein1,
                                        int protein2, GError **err)
{
  int ifeat;
  ifeat = mdt_feature_aligned_residue_add(mlib,
                                          "Residue neighborhood difference",
                                          MOD_MDTC_NONE, protein1, protein2,
                                          getbin, NULL, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_NEIGHBORS);
  }
  return ifeat;
}

int mdt_feature_average_neighborhood_difference(struct mdt_library *mlib,
                                                int protein1, int protein2,
                                                GError **err)
{
  int ifeat;
  ifeat = mdt_feature_aligned_residue_pair_add(
              mlib, "Average residue neighborhood difference", MOD_MDTC_NONE,
              protein1, protein2, FALSE, avgetbin, NULL, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_NEIGHBORS);
  }
  return ifeat;
}
