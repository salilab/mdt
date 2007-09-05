/** \file mdt_alignment.c      Functions to add alignment data to MDTs.
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
#include "mdt_property.h"

/** Get all MDT indices */
static int *mdt_indices(gboolean *outrange, const struct mod_alignment *aln,
                        int is1, int ip1, int is2, int ir1, int ir2, int ir1p,
                        int ir2p, int ia1, int ia1p,
                        const struct mdt_library *mlib, int ip2,
                        const struct mod_mdt *mdt, int ibnd1, int ibnd1p,
                        int is3, int ir3, int ir3p,
                        const struct mod_libraries *libs,
                        const struct mod_energy_data *edat,
                        struct mdt_properties *prop, GError **err)
{
  int i, *indf;
  GError *tmperr = NULL;
  *outrange = FALSE;
  indf = g_malloc(sizeof(int) * mdt->nfeat);

  for (i = 0; i < mdt->nfeat && !tmperr && *outrange == FALSE; i++) {
    const struct mod_mdt_feature *feat = &mdt->features[i];
    int ifi = feat->ifeat;
    indf[i] = my_mdt_index(ifi, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1,
                           ia1p, mlib, ip2, ibnd1, ibnd1p, is3, ir3, ir3p,
                           libs, edat, prop, &tmperr);
    if (!tmperr && (indf[i] < feat->istart || indf[i] > feat->iend)) {
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

/** Callback function to be called for every valid data point scanned for
    an MDT source. */
typedef gboolean (*mdt_scan_cb) (void *data, struct mod_mdt *mdt, int indx,
                                 GError **err);

/** Simple scan function which updates the counts in the MDT. */
static gboolean scan_update(void *data, struct mod_mdt *mdt, int indx,
                            GError **err)
{
  mdt->bin[indx] += 1.0;
  mdt->sample_size += 1.0;
  return TRUE;
}

/** Scan function which sums the MDT values. */
static gboolean scan_sum(void *data, struct mod_mdt *mdt, int indx,
                         GError **err)
{
  double *sum = (double *)data;
  *sum += mdt->bin[indx];
  return TRUE;
}

/** Update an MDT with feature data. */
static gboolean update_mdt(struct mod_mdt *mdt,
                           const struct mdt_library *mlib,
                           const struct mod_alignment *aln, int is1, int ip1,
                           int is2, int ir1, int ir2, int ir1p, int ir2p,
                           int ip2, int ia1, int ia1p, int ibnd1, int ibnd1p,
                           int is3, int ir3, int ir3p,
                           const struct mod_libraries *libs,
                           const struct mod_energy_data *edat,
                           struct mdt_properties *prop, mdt_scan_cb scanfunc,
                           void *scandata, GError **err)
{
  static const char *routine = "update_mdt";
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

  return scanfunc(scandata, mdt, imda, err);
}

/** Check to make sure that all features are single protein features */
static gboolean check_single_protein_features(const struct mod_mdt *mdt,
                                              const struct mdt_library *mlib,
                                              GError **err)
{
  int i;
  for (i = 0; i < mdt->nfeat; i++) {
    struct mod_mdt_feature *feat = &mdt->features[i];
    struct mod_mdt_libfeature *libfeat;
    libfeat = &mlib->base.features[feat->ifeat - 1];
    if (libfeat->iknown > 1) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                  "This combination of features requires a single protein "
                  "scan, but feature number %d (%d) requires multiple "
                  "proteins", i + 1, feat->ifeat);
      return FALSE;
    }
  }
  return TRUE;
}

/** Update the number of protein pairs in the MDT. */
static void update_protein_pairs(struct mod_mdt *mdt, int nseqacc,
                                 gboolean sympairs, gboolean symtriples)
{
  switch (mdt->nprotcmp) {
  case 1:
    mdt->n_protein_pairs += nseqacc;
    break;

  case 2:
    if (sympairs) {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) / 2;
    } else {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1);
    }
    break;

  case 3:
    if (symtriples) {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) * (nseqacc - 2) / 6;
    } else {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) * (nseqacc - 2);
    }
    break;
  }
}


/** Get start sequence for multiple protein features. */
static int isbeg(int is, gboolean symmetric)
{
  if (symmetric) {
    return is + 1;
  } else {
    return 0;
  }
}

/** Update MDT data for a single protein property. */
static gboolean update_single(struct mod_mdt *mdt,
                              const struct mdt_library *mlib,
                              const struct mod_alignment *aln, int is1,
                              int ip1, int ip2, int ir1, int ir1p,
                              const struct mod_libraries *libs,
                              const struct mod_energy_data *edat,
                              struct mdt_properties *prop,
                              mdt_scan_cb scanfunc, void *scandata,
                              GError **err)
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
                    ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, prop,
                    scanfunc, scandata, err)) {
      return FALSE;
    }
    break;

    /* residue properties tabulated */
  case 2:
    if (ir1 >= 0) {
      if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ip2,
                      ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, prop,
                      scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* residue rels compared */
  case 3:
    if (ir1 >= 0 && ir1p >= 0) {
      if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                      ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, prop,
                      scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;
  }
  return TRUE;
}


/** Update MDT data for a multiple protein property. */
static gboolean update_multiple(struct mod_mdt *mdt,
                                const struct mdt_library *mlib,
                                const struct mod_alignment *aln, int is1,
                                int ip1, int ip2, int ir1, int ir1p,
                                gboolean sympairs, gboolean symtriples,
                                const struct mod_libraries *libs,
                                const struct mod_energy_data *edat,
                                const gboolean acceptd[],
                                struct mdt_properties *prop,
                                mdt_scan_cb scanfunc, void *scandata,
                                GError **err)
{
  int is2, ir2, ir2p, ia1, ia1p, is3, ir3, ir3p;

  ia1 = ia1p = 0;

  /* generate all indices for the protein B: */
  for (is2 = isbeg(is1, sympairs); is2 < aln->nseq; is2++) {
    if (acceptd[is2]) {
      /* residue indices in the first and second position for protein 2 */
      if (mdt->nresfeat != 1) {
        ir2 = mod_int2_get(&aln->ialn, ip1, is2) - 1;
      }
      if (mdt->nresfeat == 3) {
        ir2p = mod_int2_get(&aln->ialn, ip2, is2) - 1;
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
                          prop, scanfunc, scandata, err)) {
            return FALSE;
          }
        }
      } else {
        /* TRIPLET OF PROTEINS:
           generate all indices for the protein C: */
        for (is3 = isbeg(is2, symtriples); is3 < aln->nseq; is3++) {
          if (acceptd[is3]) {
            if (mdt->nresfeat != 1) {
              ir3 = mod_int2_get(&aln->ialn, ip1, is3) - 1;
            }
            if (mdt->nresfeat == 3) {
              ir3p = mod_int2_get(&aln->ialn, ip2, is3) - 1;
            }
            if ((is1 != is2 && is1 != is3 && is2 != is3) || aln->nseq == 1) {
              if (!update_mdt(mdt, mlib, aln, is1, ip1, is2, ir1, ir2, ir1p,
                              ir2p, ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs,
                              edat, prop, scanfunc, scandata, err)) {
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
static gboolean genpair(struct mod_mdt *mdt, const struct mdt_library *mlib,
                        const struct mod_alignment *aln, int ip1, int ip2,
                        const struct mod_libraries *libs,
                        const struct mod_energy_data *edat,
                        const gboolean acceptd[], gboolean sympairs,
                        gboolean symtriples, struct mdt_properties *prop,
                        mdt_scan_cb scanfunc, void *scandata, GError **err)
{
  int is1, ir1, ir1p;

  /* generate all indices for protein A: */
  for (is1 = 0; is1 < aln->nseq; is1++) {
    if (acceptd[is1]) {

      /* residue index for a residue of protein A in the 1st position: */
      if (mdt->nresfeat != 1) {
        ir1 = mod_int2_get(&aln->ialn, ip1, is1) - 1;
      }
      /* residue index for a residue of protein A in the 2nd position:
         (not used if residue relationships are not compared) */
      if (mdt->nresfeat == 3 || mdt->nresfeat == 5) {
        ir1p = mod_int2_get(&aln->ialn, ip2, is1) - 1;
      }

      if (mdt->nprotcmp == 1) {
        if (!update_single
            (mdt, mlib, aln, is1, ip1, ip2, ir1, ir1p, libs, edat, prop,
             scanfunc, scandata, err)) {
          return FALSE;
        }
      } else {
        if (!update_multiple(mdt, mlib, aln, is1, ip1, ip2, ir1, ir1p, sympairs,
                             symtriples, libs, edat, acceptd, prop, scanfunc,
                             scandata, err)) {
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}


/** Return TRUE iff the two residues satisfy the residue_span_range check */
static gboolean check_sequence_separation(int ir1, int ir1p,
                                          const int rsrang[4])
{
  int nr = ir1p - ir1;
  return ((nr >= rsrang[0] && nr <= rsrang[1])
          || (nr >= rsrang[2] && nr <= rsrang[3]));
}


/** Scan all residue pairs in the first alignment sequence. */
static gboolean gen_residue_pairs(struct mod_mdt *mdt,
                                  const struct mdt_library *mlib,
                                  const struct mod_alignment *aln,
                                  const int rsrang[4],
                                  const struct mod_libraries *libs,
                                  const struct mod_energy_data *edat,
                                  const gboolean acceptd[], gboolean sympairs,
                                  gboolean symtriples,
                                  struct mdt_properties *prop,
                                  mdt_scan_cb scanfunc, void *scandata,
                                  GError **err)
{
  int ip1, ip2;

  for (ip1 = 0; ip1 < aln->naln - 1; ip1++) {

    /* only if any of the residue relationships is asymmetric, go NxN
       (mainchain H-bonds are an example) */
    if (mdt->symmetric) {
      for (ip2 = ip1 + rsrang[2]; ip2 <= MIN(aln->naln - 1, ip1 + rsrang[3]);
           ip2++) {
        if (!genpair(mdt, mlib, aln, ip1, ip2, libs, edat, acceptd, sympairs,
                     symtriples, prop, scanfunc, scandata, err)) {
          return FALSE;
        }
      }
    } else {
      for (ip2 = MAX(0, ip1 + rsrang[0]);
           ip2 <= MIN(aln->naln - 1, ip1 + rsrang[3]); ip2++) {
        if (check_sequence_separation(ip1, ip2, rsrang)) {
          if (!genpair(mdt, mlib, aln, ip1, ip2, libs, edat, acceptd, sympairs,
                       symtriples, prop, scanfunc, scandata, err)) {
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all atoms in the first alignment sequence. */
static gboolean gen_atoms(struct mod_mdt *mdt, const struct mdt_library *mlib,
                          const struct mod_alignment *aln, int is1,
                          const struct mod_libraries *libs,
                          const struct mod_energy_data *edat,
                          struct mdt_properties *prop, mdt_scan_cb scanfunc,
                          void *scandata, GError **err)
{
  int ia1, ir1, *iresatm;
  struct mod_structure *s1;

  s1 = mod_alignment_structure_get(aln, is1);

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, 1, 1, 1, ia1, 1, 1, 1,
                    1, 1, 1, libs, edat, prop, scanfunc, scandata, err)) {
      return FALSE;
    }
  }
  return TRUE;
}

/** Scan all atom pairs in the first alignment sequence. */
static gboolean gen_atom_pairs(struct mod_mdt *mdt,
                               const struct mdt_library *mlib,
                               const struct mod_alignment *aln,
                               const int rsrang[4], int is1,
                               const struct mod_libraries *libs,
                               const struct mod_energy_data *edat,
                               struct mdt_properties *prop,
                               mdt_scan_cb scanfunc, void *scandata,
                               GError **err)
{
  int ia1, ia1p, ir1, ir1p, *iresatm;
  struct mod_structure *s1;

  s1 = mod_alignment_structure_get(aln, is1);

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ia1p = ia1 + 1; ia1p < s1->cd.natm; ia1p++) {
      ir1p = iresatm[ia1p] - 1;
      if (check_sequence_separation(ir1, ir1p, rsrang)) {
        if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1,
                        ia1p, 1, 1, 1, 1, 1, libs, edat, prop, scanfunc,
                        scandata, err)) {
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}


/** Scan all bonds, angles or dihedrals in the first alignment sequence. */
static gboolean gen_bonds(struct mod_mdt *mdt, const struct mdt_library *mlib,
                          const struct mod_alignment *aln, int is1, int npnt,
                          const struct mod_libraries *libs,
                          const struct mod_energy_data *edat,
                          struct mdt_properties *prop, mdt_scan_cb scanfunc,
                          void *scandata, GError **err)
{
  const struct mdt_bond_list *bonds;
  struct mod_structure *struc;
  int ibnd1, is2;

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  struc = mod_alignment_structure_get(aln, is1);
  is2 = is1;
  bonds = property_bonds(aln, is1, prop, mlib, npnt, libs);
  for (ibnd1 = 0; ibnd1 < bonds->nbonds; ibnd1++) {
    if (!update_mdt(mdt, mlib, aln, is1, 1, is2, 1, 1, 1, 1, 1, 1, 1, ibnd1,
                    1, 1, 1, 1, libs, edat, prop, scanfunc, scandata, err)) {
      return FALSE;
    }
  }
  return TRUE;
}


/** Scan all atom tuples in the first alignment sequence. */
static gboolean gen_atom_tuples(struct mod_mdt *mdt,
                                const struct mdt_library *mlib,
                                const struct mod_alignment *aln, int is1,
                                const struct mod_libraries *libs,
                                const struct mod_energy_data *edat,
                                struct mdt_properties *prop,
                                mdt_scan_cb scanfunc, void *scandata,
                                GError **err)
{
  int ia1, ir1, ibnd1, ibnd1p, ia1p, ir1p, *iresatm;
  struct mod_structure *s1;
  const struct mdt_tuple_list *tup;

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  s1 = mod_alignment_structure_get(aln, is1);
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  tup = property_tuples(aln, is1, prop, mlib, libs);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ibnd1 = 0; ibnd1 < tup[ia1].ntuples; ibnd1++) {
      /* Just in case you use a single atom feature at position 2 in
         protein A: */
      ia1p = ia1;
      ibnd1p = ibnd1;
      ir1p = ir1;
      if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1, ia1p,
                      ibnd1, ibnd1p, 1, 1, 1, libs, edat, prop, scanfunc,
                      scandata, err)) {
        return FALSE;
      }
    }
  }
  return TRUE;
}


/** Scan all atom tuple pairs in the first alignment sequence. */
static gboolean gen_atom_tuple_pairs(struct mod_mdt *mdt,
                                     const struct mdt_library *mlib,
                                     const struct mod_alignment *aln,
                                     const int rsrang[4], int is1,
                                     const struct mod_libraries *libs,
                                     const struct mod_energy_data *edat,
                                     struct mdt_properties *prop,
                                     mdt_scan_cb scanfunc, void *scandata,
                                     GError **err)
{
  int ia1, ir1, ibnd1, ibnd1p, ia1p, ir1p, nr, *iresatm;
  struct mod_structure *s1;
  const struct mdt_tuple_list *tup;

  s1 = mod_alignment_structure_get(aln, is1);
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  tup = property_tuples(aln, is1, prop, mlib, libs);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ibnd1 = 0; ibnd1 < tup[ia1].ntuples; ibnd1++) {
      for (ia1p = 0; ia1p < s1->cd.natm; ia1p++) {
        ir1p = iresatm[ia1p] - 1;

        /* the same conditions on sequence separation as for residue pairs */
        nr = ir1p - ir1;
        if (ia1 != ia1p && check_sequence_separation(ir1, ir1p, rsrang)) {
          for (ibnd1p = 0; ibnd1p < tup[ia1p].ntuples; ibnd1p++) {
            if (!update_mdt(mdt, mlib, aln, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1,
                            ia1p, ibnd1, ibnd1p, 1, 1, 1, libs, edat, prop,
                            scanfunc, scandata, err)) {
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
static gboolean mdt_source_scan(struct mod_mdt *mdt,
                                const struct mdt_library *mlib,
                                const struct mod_alignment *aln,
                                const int rsrang[4],
                                const struct mod_libraries *libs,
                                const struct mod_energy_data *edat,
                                const gboolean acceptd[], int nseqacc,
                                gboolean sympairs, gboolean symtriples,
                                struct mdt_properties *prop,
                                mdt_scan_cb scanfunc, void *scandata,
                                GError **err)
{
  int ip1;

  update_protein_pairs(mdt, nseqacc, sympairs, symtriples);

  switch (mdt->nresfeat) {
    /* Whole proteins */
  case 1:
    if (!genpair(mdt, mlib, aln, 1, 1, libs, edat, acceptd, sympairs,
                 symtriples, prop, scanfunc, scandata, err)) {
      return FALSE;
    }
    break;

    /* Single residues or selected (one per residue) atoms */
  case 2:
  case 4:
    for (ip1 = 0; ip1 < aln->naln; ip1++) {
      if (!genpair(mdt, mlib, aln, ip1, ip1, libs, edat, acceptd, sympairs,
                   symtriples, prop, scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* intra-molecular residue or selected (one per residue) atom pairs */
  case 3:
  case 5:
    if (!gen_residue_pairs(mdt, mlib, aln, rsrang, libs, edat, acceptd,
                           sympairs, symtriples, prop, scanfunc, scandata,
                           err)) {
      return FALSE;
    }
    break;

    /* Single protein, all atoms; using only the first protein in
       an alignment! */
  case 6:
    if (acceptd[0]) {
      if (!gen_atoms(mdt, mlib, aln, 0, libs, edat, prop, scanfunc, scandata,
                     err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom pairs; using only the first protein in
       an alignment! */
  case 7:
    if (acceptd[0]) {
      if (!gen_atom_pairs(mdt, mlib, aln, rsrang, 0, libs, edat, prop,
                          scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom tuples; using only the first protein in
       an alignment! */
  case 8:
    if (acceptd[0]) {
      if (!gen_atom_tuples(mdt, mlib, aln, 0, libs, edat, prop, scanfunc,
                           scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom tuple pairs; using only the first protein in
       an alignment! */
  case 9:
    if (acceptd[0]) {
      if (!gen_atom_tuple_pairs(mdt, mlib, aln, rsrang, 0, libs, edat, prop,
                                scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all bonds: */
  case 10:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, aln, 0, MDT_BOND_TYPE_BOND, libs, edat, prop,
                     scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all angles: */
  case 11:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, aln, 0, MDT_BOND_TYPE_ANGLE, libs, edat, prop,
                     scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all dihedrals: */
  case 12:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, aln, 0, MDT_BOND_TYPE_DIHEDRAL, libs, edat,
                     prop, scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;
  }
  return TRUE;
}


/** A source of data for an MDT (generally an alignment) */
struct mdt_source {
  int nseqacc;
  struct mdt_properties *prop;
  struct mod_alignment *aln;
  gboolean *acceptd;
};


/** Prepare a source alignment to add data to an MDT. Returns a source pointer
    (to be later freed with mdt_alignment_close()), or NULL on error. */
struct mdt_source *mdt_alignment_open(struct mod_mdt *mdt,
                                      const struct mdt_library *mlib,
                                      struct mod_alignment *aln, float distngh,
                                      gboolean sdchngh, int surftyp,
                                      int iacc1typ, struct mod_io_data *io,
                                      struct mod_libraries *libs, GError **err)
{
  int nseqacc, ierr;

  mod_mdt_getdata(mdt, &nseqacc, aln, distngh, sdchngh, surftyp, iacc1typ, io,
                  libs, &ierr);
  if (ierr) {
    handle_modeller_error(err);
    return NULL;
  }

  mod_lognote("Pre-calculating");
  mod_mdt_precalc(mdt, &mlib->base, aln, libs, &ierr);
  if (ierr) {
    handle_modeller_error(err);
    return NULL;
  } else {
    struct mdt_source *source = g_malloc(sizeof(struct mdt_source));
    source->nseqacc = nseqacc;
    source->aln = aln;
    source->prop = mdt_properties_new(aln);
    source->acceptd = g_malloc(sizeof(gboolean) * aln->nseq);

    if (mdt->readin[0]) {
      int i;
      for (i = 0; i < aln->nseq; i++) {
        struct mod_structure *s1 = mod_alignment_structure_get(aln, i);
        source->acceptd[i] = s1->accepts;
      }
    } else {
      int i;
      for (i = 0; i < aln->nseq; i++) {
        source->acceptd[i] = TRUE;
      }
    }
    return source;
  }
}


/** Close a source alignment previously opened with mdt_alignment_open(). */
void mdt_alignment_close(struct mdt_source *source)
{
  mdt_properties_free(source->prop, source->aln);
  g_free(source->acceptd);
  g_free(source);
}


/** Scan all data points in the source, and return the sum. */
double mdt_source_sum(struct mdt_source *source, struct mod_mdt *mdt,
                      const struct mdt_library *mlib,
                      const int residue_span_range[4],
                      const struct mod_libraries *libs,
                      const struct mod_energy_data *edat, gboolean sympairs,
                      gboolean symtriples, GError **err)
{
  double sum = 0.;
  mdt_source_scan(mdt, mlib, source->aln, residue_span_range, libs,
                  edat, source->acceptd, source->nseqacc, sympairs, symtriples,
                  source->prop, scan_sum, &sum, err);
  return sum;
}


/** Return the bin index (starting at 1) of a single MDT feature, at the
    given position in the source alignment. On failure, 0 is returned. */
int mdt_alignment_index(struct mdt_source *source, int ifeat, int is1, int ip1,
                        int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                        int ia1p, int ip2, int ibnd1, int ibnd1p, int is3,
                        int ir3, int ir3p, const struct mdt_library *mlib,
                        const struct mod_libraries *libs,
                        struct mod_energy_data *edat, GError **err)
{
  int indf;
  indf = my_mdt_index(ifeat, source->aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                      ia1, ia1p, mlib, ip2, ibnd1, ibnd1p, is3, ir3,
                      ir3p, libs, edat, source->prop, err);
  return indf;
}


/** Add data from an alignment to an MDT. Return TRUE on success. */
gboolean mdt_add_alignment(struct mod_mdt *mdt,
                           const struct mdt_library *mlib,
                           struct mod_alignment *aln, float distngh,
                           gboolean sdchngh, int surftyp, int iacc1typ,
                           const int residue_span_range[4], gboolean sympairs,
                           gboolean symtriples, struct mod_io_data *io,
                           struct mod_energy_data *edat,
                           struct mod_libraries *libs, GError **err)
{
  struct mdt_source *source;

  mod_lognote("Calculating and checking other data: %d", aln->nseq);

  source = mdt_alignment_open(mdt, mlib, aln, distngh, sdchngh, surftyp,
                              iacc1typ, io, libs, err);
  if (source) {
    gboolean ret;

    mdt->nalns++;
    mdt->n_proteins += source->nseqacc;

    mod_lognote("Updating the statistics array:");
    ret = mdt_source_scan(mdt, mlib, aln, residue_span_range, libs, edat,
                          source->acceptd, source->nseqacc, sympairs,
                          symtriples, source->prop, scan_update, NULL, err);

    mdt_alignment_close(source);
    return ret;
  } else {
    return FALSE;
  }
}
