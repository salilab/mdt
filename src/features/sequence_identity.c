/** \file sequence_identity.c  Sequence identity feature.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

/** Return the number of identical aligned residues in the two sequences */
static int num_equiv(const struct mod_alignment *aln, int protein1,
                     int protein2, const struct mod_sequence *s1,
                     const struct mod_sequence *s2)
{
  int i, neqv = 0;

  for (i = 0; i < aln->naln; ++i) {
    int ires1 = mod_int2_get(&aln->ialn, i, protein1) - 1;
    int ires2 = mod_int2_get(&aln->ialn, i, protein2) - 1;
    if (ires1 >= 0 && ires2 >= 0) {
      int restyp1 = mod_int1_get(&s1->irestyp, ires1);
      int restyp2 = mod_int1_get(&s2->irestyp, ires2);
      if (restyp1 == restyp2) {
        neqv++;
      }
    }
  }
  return neqv;
}

static int getbin(const struct mod_alignment *aln, int protein1, int protein2,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_sequence *s1 = mod_alignment_sequence_get(aln, protein1);
  struct mod_sequence *s2 = mod_alignment_sequence_get(aln, protein2);
  int nres;
  float f;
  f = num_equiv(aln, protein1, protein2, s1, s2);
  nres = MIN(s1->nres, s2->nres);
  if (nres > 0) {
    f /= nres;
  }
  return feat_to_bin(f, feat);
}

int mdt_feature_sequence_identity(struct mdt_library *mlib, int protein1,
                                  int protein2, GError **err)
{
  return mdt_feature_protein_pair_add(mlib,
                                      "Overall fractional sequence identity",
                                      MOD_MDTC_NONE, protein1, protein2, getbin,
                                      NULL, NULL, err);
}
