/** \file mdt_add_alignment.c  Functions to add alignment data to MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include <assert.h>
#include <glib.h>

#include "modeller.h"
#include "mod_error.h"
#include "mdt.h"
#include "util.h"
#include "mdt_index.h"

/** Get all MDT indices */
static int *mdt_indices(gboolean *outrange, const struct alignment *aln,
                        int is1, int ip1, int is2, int ir1, int ir2, int ir1p,
                        int ir2p, int ia1, int ia1p,
                        const struct mdt_library *mlib, int ip2,
                        const struct mdt_type *mdt, int ibnd1, int ibnd1p,
                        int is3, int ir3, int ir3p,
                        const struct libraries *libs,
                        const struct energy_data *edat,
                        struct mdt_properties *prop, GError **err)
{
  int i, *indf, *ifeat, *istart, *iend;
  GError *tmperr = NULL;
  *outrange = FALSE;
  indf = g_malloc(sizeof(int) * mdt->nfeat);

  ifeat = f_int1_pt(&mdt->ifeat);
  istart = f_int1_pt(&mdt->istart);
  iend = f_int1_pt(&mdt->iend);
  for (i = 0; i < mdt->nfeat && !tmperr && *outrange == FALSE; i++) {
    int ifi = ifeat[i];
    indf[i] = my_mdt_index(ifi, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1,
                           ia1p, mlib, ip2, ibnd1, ibnd1p, is3, ir3, ir3p,
                           libs, edat, prop, &tmperr);
    if (!tmperr && (indf[i] < istart[i] || indf[i] > iend[i])) {
      *outrange = TRUE;
    }
  }

  if (tmperr) {
    g_propagate_error(err, tmperr);
    return NULL;
  } else {
    return indf;
  }
}


/** Update an MDT with feature data. */
static gboolean update_mdt(struct mdt_type *mdt,
                           const struct mdt_library *mlib,
                           const struct alignment *aln, int is1, int ip1,
                           int is2, int ir1, int ir2, int ir1p, int ir2p,
                           int ip2, int ia1, int ia1p, int ibnd1, int ibnd1p,
                           int is3, int ir3, int ir3p,
                           const struct libraries *libs,
                           const struct energy_data *edat,
                           struct mdt_properties *prop, GError **err)
{
  static const char *routine = "update_mdt";
  double *bin;
  gboolean outrange;
  int imda, *indf;

  /* obtain the indices for the feature values in this routine call: */
  indf = mdt_indices(&outrange, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1,
                     ia1p, mlib, ip2, mdt, ibnd1, ibnd1p, is3, ir3, ir3p, libs,
                     edat, prop, err);
  if (!indf) {
    return FALSE;
  }

  /* Ignore if any of the indices properly out of range: */
  if (outrange) {
    g_free(indf);
    return TRUE;
  }

  /* obtain the element index for the mdt vector: */
  imda = indmdt(indf, mdt);
  g_free(indf);
  if (imda < 0 || imda >= mdt->nelems) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                "%s: MDT index is out of range: %d %d", routine, imda,
                mdt->nelems);
    return FALSE;
  }

  bin = f_double1_pt(&mdt->bin);
  bin[imda] += 1.0;
  mdt->sample_size += 1.0;
  return TRUE;
}

/** Update the number of protein pairs in the MDT. */
static void update_protein_pairs(struct mdt_type *mdt, int nseqacc, int pairs,
                                 int triples)
{
  switch (mdt->nprotcmp) {
  case 1:
    mdt->n_protein_pairs += nseqacc;
    break;

  case 2:
    if (pairs == 1) {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1);
    } else {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) / 2;
    }
    break;

  case 3:
    if (triples == 1) {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) * (nseqacc - 2);
    } else {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) * (nseqacc - 2) / 6;
    }
    break;
  }
}


/** Get start sequence for multiple protein features. */
static int isbeg(int is, int nseq, int iseqbeg)
{
  switch (iseqbeg) {
  default:
    return 0;
  case 2:
    return is + 1;
  case 3:
    return nseq - 1;
  }
}

/** Update MDT data for a single protein property. */
static gboolean update_single(struct mdt_type *mdt,
                              const struct mdt_library *mlib,
                              const struct alignment *aln, int is1, int ip1,
                              int ip2, int ir1, int ir1p,
                              const struct libraries *libs,
                              const struct energy_data *edat,
                              struct mdt_properties *prop, GError **err)
{
  int is2, ir2, ir2p, ia1, ia1p, is3, ir3, ir3p;

  is2 = is3 = is1;
  ir2 = ir3 = ir1;
  ir2p = ir3p = ir1p;
  ia1 = ia1p = 0;

  switch (mdt->nresfeat) {
    /* whole protein properties tabulated: */
  case 1:
    if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ip2,
                    ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, prop, err)) {
      return FALSE;
    }
    break;

    /* residue properties tabulated */
  case 2:
    if (ir1 >= 0) {
      if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ip2,
                      ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, prop,
                      err)) {
        return FALSE;
      }
    }
    break;

    /* residue rels compared */
  case 3:
    if (ir1 >= 0 && ir1p >= 0) {
      if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                      ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, prop,
                      err)) {
        return FALSE;
      }
    }
    break;
  }
  return TRUE;
}


/** Update MDT data for a multiple protein property. */
static gboolean update_multiple(struct mdt_type *mdt,
                                const struct mdt_library *mlib,
                                const struct alignment *aln, int is1, int ip1,
                                int ip2, int ir1, int ir1p, int pairs,
                                int triples, const struct libraries *libs,
                                const struct energy_data *edat,
                                const gboolean acceptd[],
                                struct mdt_properties *prop, GError **err)
{
  int is2, ir2, ir2p, ia1, ia1p, is3, ir3, ir3p;

  ia1 = ia1p = 0;

  /* generate all indices for the protein B: */
  for (is2 = isbeg(is1, aln->nseq, pairs); is2 < aln->nseq; is2++) {
    if (acceptd[is2]) {
      /* residue indices in the first and second position for protein 2 */
      if (mdt->nresfeat != 1) {
        ir2 = f_int2_get(&aln->ialn, ip1, is2) - 1;
      }
      if (mdt->nresfeat == 3) {
        ir2p = f_int2_get(&aln->ialn, ip2, is2) - 1;
      }
      is3 = is2;
      ir3 = ir2;
      ir3p = ir2p;

      /* do you have to generate a pair or a triplet of proteins */
      if (mdt->nprotcmp == 2) {
        /* PAIR OF PROTEINS:
           ignore self-comparisons:
           (the second condition only applies in the GETCSR mode to
           allow the prediction from the sequence of the unknown alone
           even if the MDT table was read in that requires the known) */
        if (is1 != is2 || aln->nseq == 1) {
          if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                          ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat,
                          prop, err)) {
            return FALSE;
          }
        }
      } else {
        /* TRIPLET OF PROTEINS:
           generate all indices for the protein C: */
        for (is3 = isbeg(is2, aln->nseq, triples); is3 < aln->nseq; is3++) {
          if (acceptd[is3]) {
            if (mdt->nresfeat != 1) {
              ir3 = f_int2_get(&aln->ialn, ip1, is3) - 1;
            }
            if (mdt->nresfeat == 3) {
              ir3p = f_int2_get(&aln->ialn, ip2, is3) - 1;
            }
            if ((is1 != is2 && is1 != is3 && is2 != is3) || aln->nseq == 1) {
              if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p,
                              ir2p, ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs,
                              edat, prop, err)) {
                return FALSE;
              }
            }
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all proteins or protein pairs in the alignment. */
static gboolean genpair(struct mdt_type *mdt, const struct mdt_library *mlib,
                        const struct alignment *aln, int ip1, int ip2,
                        const struct libraries *libs,
                        const struct energy_data *edat,
                        const gboolean acceptd[], int pairs, int triples,
                        struct mdt_properties *prop, GError **err)
{
  int is1, ir1, ir1p;

  /* generate all indices for protein A: */
  for (is1 = 0; is1 < aln->nseq; is1++) {
    if (acceptd[is1]) {

      /* residue index for a residue of protein A in the 1st position: */
      if (mdt->nresfeat != 1) {
        ir1 = f_int2_get(&aln->ialn, ip1, is1) - 1;
      }
      /* residue index for a residue of protein A in the 2nd position:
         (not used if residue relationships are not compared) */
      if (mdt->nresfeat == 3 || mdt->nresfeat == 5) {
        ir1p = f_int2_get(&aln->ialn, ip2, is1) - 1;
      }

      if (mdt->nprotcmp == 1) {
        if (!update_single
            (mdt, mlib, aln, is1, ip1, ip2, ir1, ir1p, libs, edat, prop,
             err)) {
          return FALSE;
        }
      } else {
        if (!update_multiple(mdt, mlib, aln, is1, ip1, ip2, ir1, ir1p, pairs,
                             triples, libs, edat, acceptd, prop, err)) {
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}



/** Scan all residue pairs in the first alignment sequence. */
static gboolean gen_residue_pairs(struct mdt_type *mdt,
                                  const struct mdt_library *mlib,
                                  const struct alignment *aln,
                                  const int rsrang[4],
                                  const struct libraries *libs,
                                  const struct energy_data *edat,
                                  const gboolean acceptd[], int pairs,
                                  int triples, struct mdt_properties *prop,
                                  GError **err)
{
  int ip1, ip2;

  for (ip1 = 0; ip1 < aln->naln - 1; ip1++) {

    /* only if any of the residue relationships is asymmetric, go NxN
       (mainchain H-bonds are an example) */
    if (mdt->symmetric) {
      for (ip2 = ip1 + rsrang[2]; ip2 <= MIN(aln->naln - 1, ip1 + rsrang[3]);
           ip2++) {
        if (!genpair(mdt, mlib, aln, ip1, ip2, libs, edat, acceptd, pairs,
                     triples, prop, err)) {
          return FALSE;
        }
      }
    } else {
      for (ip2 = MAX(0, ip1 - rsrang[3]);
           ip2 <= MIN(aln->naln - 1, ip1 + rsrang[3]); ip2++) {
        if (abs(ip1 - ip2) >= rsrang[1]) {
          if (!genpair(mdt, mlib, aln, ip1, ip2, libs, edat, acceptd, pairs,
                       triples, prop, err)) {
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all atoms in the first alignment sequence. */
static gboolean gen_atoms(struct mdt_type *mdt, const struct mdt_library *mlib,
                          const struct alignment *aln, int is1,
                          const struct libraries *libs,
                          const struct energy_data *edat,
                          struct mdt_properties *prop, GError **err)
{
  int ia1, ir1, *iresatm;
  struct structure *s1;

  s1 = alignment_structure_get(aln, is1);

  iresatm = f_int1_pt(&s1->cd.iresatm);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, 1, 1, 1, ia1, 1, 1, 1,
                    1, 1, 1, libs, edat, prop, err)) {
      return FALSE;
    }
  }
  return TRUE;
}

/** Scan all atom pairs in the first alignment sequence. */
static gboolean gen_atom_pairs(struct mdt_type *mdt,
                               const struct mdt_library *mlib,
                               const struct alignment *aln, int is1,
                               const struct libraries *libs,
                               const struct energy_data *edat,
                               struct mdt_properties *prop, GError **err)
{
  int ia1, ia1p, ir1, ir1p, *iresatm;
  struct structure *s1;

  s1 = alignment_structure_get(aln, is1);

  iresatm = f_int1_pt(&s1->cd.iresatm);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ia1p = ia1 + 1; ia1p < s1->cd.natm; ia1p++) {
      ir1p = iresatm[ia1p] - 1;
      if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1,
                      ia1p, 1, 1, 1, 1, 1, libs, edat, prop, err)) {
        return FALSE;
      }
    }
  }
  return TRUE;
}


/** Scan all bonds, angles or dihedrals in the first alignment sequence. */
static gboolean gen_bonds(struct mdt_type *mdt, const struct mdt_library *mlib,
                          const struct alignment *aln, int is1, int npnt,
                          const struct libraries *libs,
                          const struct energy_data *edat,
                          struct mdt_properties *prop, GError **err)
{
  const struct mdt_bond_list *bonds;
  struct structure *struc;
  int ibnd1, is2;

  struc = alignment_structure_get(aln, is1);
  is2 = is1;
  bonds = property_bonds(aln, is1, prop, mlib, npnt, libs);
  for (ibnd1 = 0; ibnd1 < bonds->nbonds; ibnd1++) {
    if (!update_mdt(mdt, mlib, aln, is1, 1, is2, 1, 1, 1, 1, 1, 1, 1, ibnd1,
                    1, 1, 1, 1, libs, edat, prop, err)) {
      return FALSE;
    }
  }
  return TRUE;
}


/** Scan all atom triplets in the first alignment sequence. */
static gboolean gen_atom_triplets(struct mdt_type *mdt,
                                  const struct mdt_library *mlib,
                                  const struct alignment *aln, int is1,
                                  const struct libraries *libs,
                                  const struct energy_data *edat,
                                  struct mdt_properties *prop, GError **err)
{
  int ia1, ir1, ibnd1, ibnd1p, ia1p, ir1p, *iresatm;
  struct structure *s1;
  const struct mdt_triplet_list *trp;

  s1 = alignment_structure_get(aln, is1);
  iresatm = f_int1_pt(&s1->cd.iresatm);
  trp = property_triplets(aln, is1, prop, mlib, libs);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ibnd1 = 0; ibnd1 < trp[ia1].ntriplets; ibnd1++) {
      /* Just in case you use a single atom feature at position 2 in
         protein A: */
      ia1p = ia1;
      ibnd1p = ibnd1;
      ir1p = ir1;
      if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1, ia1p,
                      ibnd1, ibnd1p, 1, 1, 1, libs, edat, prop, err)) {
        return FALSE;
      }
    }
  }
  return TRUE;
}


/** Scan all atom triplet pairs in the first alignment sequence. */
static gboolean gen_atom_triplet_pairs(struct mdt_type *mdt,
                                       const struct mdt_library *mlib,
                                       const struct alignment *aln,
                                       const int rsrang[4], int is1,
                                       const struct libraries *libs,
                                       const struct energy_data *edat,
                                       struct mdt_properties *prop,
                                       GError **err)
{
  int ia1, ir1, ibnd1, ibnd1p, ia1p, ir1p, nr, *iresatm;
  struct structure *s1;
  const struct mdt_triplet_list *trp;

  s1 = alignment_structure_get(aln, is1);
  iresatm = f_int1_pt(&s1->cd.iresatm);
  trp = property_triplets(aln, is1, prop, mlib, libs);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ibnd1 = 0; ibnd1 < trp[ia1].ntriplets; ibnd1++) {
      for (ia1p = 0; ia1p < s1->cd.natm; ia1p++) {
        ir1p = iresatm[ia1p] - 1;

        /* the same conditions on sequence separation as for residue pairs */
        nr = ir1p - ir1;
        if (ia1 != ia1p && ((nr >= rsrang[0] && nr <= rsrang[1])
                            || (nr >= rsrang[2] && nr <= rsrang[3]))) {
          for (ibnd1p = 0; ibnd1p < trp[ia1p].ntriplets; ibnd1p++) {
            if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1,
                            ia1p, ibnd1, ibnd1p, 1, 1, 1, libs, edat, prop,
                            err)) {
              return FALSE;
            }
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all alignment positions or all alignment position pairs in the
    current alignment. If whole protein features only occur in the current
    MDT, then no positions are scanned. */
static gboolean update_stats(struct mdt_type *mdt,
                             const struct mdt_library *mlib,
                             const struct alignment *aln, const int rsrang[4],
                             const struct libraries *libs,
                             const struct energy_data *edat,
                             const gboolean acceptd[], int nseqacc, int pairs,
                             int triples, struct mdt_properties *prop,
                             GError **err)
{
  int ip1;

  update_protein_pairs(mdt, nseqacc, pairs, triples);

  switch (mdt->nresfeat) {
    /* Whole proteins */
  case 1:
    if (!genpair(mdt, mlib, aln, 1, 1, libs, edat, acceptd, pairs, triples,
                 prop, err)) {
      return FALSE;
    }
    break;

    /* Single residues or selected (one per residue) atoms */
  case 2:
  case 4:
    for (ip1 = 0; ip1 < aln->naln; ip1++) {
      if (!genpair(mdt, mlib, aln, ip1, ip1, libs, edat, acceptd, pairs,
                   triples, prop, err)) {
        return FALSE;
      }
    }
    break;

    /* intra-molecular residue or selected (one per residue) atom pairs */
  case 3:
  case 5:
    if (!gen_residue_pairs(mdt, mlib, aln, rsrang, libs, edat, acceptd, pairs,
                           triples, prop, err)) {
      return FALSE;
    }
    break;

    /* Single protein, all atoms; using only the first protein in
       an alignment! */
  case 6:
    if (acceptd[0]) {
      if (!gen_atoms(mdt, mlib, aln, 0, libs, edat, prop, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom pairs; using only the first protein in
       an alignment! */
  case 7:
    if (acceptd[0]) {
      if (!gen_atom_pairs(mdt, mlib, aln, 0, libs, edat, prop, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom triplets; using only the first protein in
       an alignment! */
  case 8:
    if (acceptd[0]) {
      if (!gen_atom_triplets(mdt, mlib, aln, 0, libs, edat, prop, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom triplet pairs; using only the first protein in
       an alignment! */
  case 9:
    if (acceptd[0]) {
      if (!gen_atom_triplet_pairs(mdt, mlib, aln, rsrang, 0, libs, edat, prop,
                                  err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all bonds: */
  case 10:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, aln, 0, MDT_BOND_TYPE_BOND, libs, edat, prop,
                     err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all angles: */
  case 11:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, aln, 0, MDT_BOND_TYPE_ANGLE, libs, edat, prop,
                     err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all dihedrals: */
  case 12:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, aln, 0, MDT_BOND_TYPE_DIHEDRAL, libs, edat,
                     prop, err)) {
        return FALSE;
      }
    }
    break;
  }
  return TRUE;
}


/** Add data from an alignment to an MDT. Return TRUE on success. */
gboolean mdt_add_alignment(struct mdt_type *mdt,
                           const struct mdt_library *mlib,
                           struct alignment *aln, float distngh,
                           gboolean sdchngh, int surftyp, int iacc1typ,
                           const int residue_span_range[4], int pairs,
                           int triples, struct io_data *io,
                           struct energy_data *edat, struct libraries *libs,
                           GError **err)
{
  int nseqacc, ierr;
  gboolean ret, *acceptd;

  modlognote("Calculating and checking other data: %d", aln->nseq);

  mdt_getdata(mdt, &nseqacc, aln, distngh, sdchngh, surftyp, iacc1typ, io,
              libs, &ierr);
  if (ierr) {
    handle_modeller_error(err);
    return FALSE;
  }

  acceptd = g_malloc(sizeof(gboolean) * aln->nseq);
  if (mdt->readin[0]) {
    int i;
    for (i = 0; i < aln->nseq; i++) {
      struct structure *s1 = alignment_structure_get(aln, i);
      acceptd[i] = s1->accepts;
    }
  } else {
    int i;
    for (i = 0; i < aln->nseq; i++) {
      acceptd[i] = TRUE;
    }
  }

  mdt->nalns++;
  mdt->n_proteins += nseqacc;

  modlognote("Pre-calculating");
  mdt_precalc(mdt, &mlib->base, aln, libs, &ierr);
  if (ierr) {
    handle_modeller_error(err);
    ret = FALSE;
  } else {
    struct mdt_properties *prop = mdt_properties_new(aln);
    modlognote("Updating the statistics array:");
    ret = update_stats(mdt, mlib, aln, residue_span_range, libs, edat, acceptd,
                       nseqacc, pairs, triples, prop, err);
    mdt_properties_free(prop, aln);
  }
  g_free(acceptd);
  return ret;
}
