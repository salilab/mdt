/** \file mdt_index.c  Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <string.h>
#include <math.h>
#include "modeller.h"
#include "util.h"
#include "mdt_index.h"
#include "mdt_hydrogen_bonds.h"

/** Make a new mdt_properties structure */
struct mdt_properties *mdt_properties_new(const struct alignment *aln)
{
  struct mdt_properties *prop;
  int i;
  prop = g_malloc(sizeof(struct mdt_properties) * aln->naln);
  for (i = 0; i < aln->naln; i++) {
    prop[i].hb_iatta = NULL;
    prop[i].hbpot = NULL;
    prop[i].iatta = NULL;
    prop[i].iatmacc = NULL;
    prop[i].ifatmacc = NULL;
  }
  return prop;
}

/** Free an mdt_properties structure */
void mdt_properties_free(struct mdt_properties *prop,
                         const struct alignment *aln)
{
  int i;
  for (i = 0; i < aln->naln; i++) {
    g_free(prop[i].hb_iatta);
    g_free(prop[i].hbpot);
    g_free(prop[i].iatta);
    g_free(prop[i].iatmacc);
    g_free(prop[i].ifatmacc);
  }
  g_free(prop);
}

static int ires_get(int ires, int nres, int igaptyp, const int *irestyp,
                    int ndimen)
{
  if (ires < 1 || ires > nres) {
    return igaptyp;
  } else if (irestyp[ires-1] >= 21 || irestyp[ires-1] <= 0) {
    modlogwarning("irestab", "ires, irestp: %d %d", ires, irestyp[ires-1]);
    return ndimen;
  } else {
    return irestyp[ires-1];
  }
}

static int irestab(const struct f_int2_array *ialn, int naln, int iseq,
                   const int *irestyp, int nres, int ip, int delta,
                   gboolean delta_ali, int ndimen, int igaptyp)
{
  if (delta_ali) {
    int ipos = ip + delta;
    if (ipos < 1 || ipos > naln) {
      return igaptyp;
    } else {
      int ires = f_int2_get(ialn, ipos-1, iseq-1);
      return ires_get(ires, nres, igaptyp, irestyp, ndimen);
    }
  } else {
    int ires = f_int2_get(ialn, ip-1, iseq-1);
    if (ires < 1 || ires > nres) {
      return igaptyp;
    } else {
      return ires_get(ires + delta, nres, igaptyp, irestyp, ndimen);
    }
  }
}

static int itable(const int *itab, int nr, int ir, int ndim)
{
  if (ir >= 1 && ir <= nr && itab[ir-1] >= 1 && itab[ir-1] <= ndim) {
    return itab[ir-1];
  } else {
    return ndim;
  }
}

static int iclsbin(float x, const struct mdt_library *mlib, int ifi, int nrang)
{
  int i;
  const struct mdt_feature *feat = &mlib->base.features[ifi-1];
  const struct mdt_bin *bin = feat->bins;
  for (i = 0; i < nrang; i++, bin++) {
    if (x >= bin->rang1 && x <= bin->rang2) {
      return i + 1;
    }
  }
  bin = &feat->bins[0];
  modlogwarning("iclsbin", "Undefined value; X,x1,x2,n,bin: %f %f %f %d",
                x, bin->rang1, bin->rang2, nrang + 1);
  return nrang + 1;
}

static void alliclsbin(int nvec, const float *x, int *ix,
                       const struct mdt_library *mlib, int ifi, int nrang)
{
  int i;
  for (i = 0; i < nvec; i++) {
    ix[i] = iclsbin(x[i], mlib, ifi, nrang);
  }
}

static int iatmcls(int irestyp, const char *atmnam, 
                   const struct mdt_atom_class_list *atclass,
                   const struct libraries *libs)
{
  gboolean allres, allatm;
  char *resnam;
  int iclass;

  allres = (irestyp == 0);
  allatm = (strcmp(atmnam, "*") == 0);

  if (allres) {
    resnam = g_strdup("*");
  } else {
    resnam = residue_name_from_type(irestyp, libs);
  }

  for (iclass = 0; iclass < atclass->nclass; iclass++) {
    int i;
    const struct mdt_atom_class *atc = &atclass->classes[iclass];
    for (i = 0; i < atc->ntypes; i++) {
      const struct mdt_atom_type *att = &atc->types[i];
      if ((allres || strcmp(att->names[0], "*") == 0
          || strcmp(att->names[0], resnam) == 0)
          && (allatm || strcmp(att->names[1], "*") == 0
              || strcmp(att->names[1], atmnam) == 0)) {
        g_free(resnam);
        return iclass + 1;
      }
    }
  }

  if (!residue_is_hetatm(irestyp, libs) && !mod_atom_is_hydrogen(atmnam)) {
    modlogwarning("iatmcls", "Model atom not classified: %s:%s", resnam,
                  atmnam);
  }
  g_free(resnam);
  return 0;
}

static void atmcls_special(struct structure *struc, const struct sequence *seq,
                           int iatta[],
                           const struct mdt_atom_class_list *atclass,
                           const struct mdt_library *mlib,
                           const struct libraries *libs)
{
  int i, *irestyp, *iresatm;
  iresatm = f_int1_pt(&struc->cd.iresatm);
  irestyp = f_int1_pt(&seq->irestyp);
  for (i = 0; i < struc->cd.natm; i++) {
    int irest = irestyp[iresatm[i] - 1];
    char *atmnam = get_coord_atmnam(&struc->cd, i);
    iatta[i] = iatmcls(irest, atmnam, atclass, libs);
    g_free(atmnam);
  }

  if (mlib->base.special_atoms) {
    int *iatmr1 = f_int1_pt(&struc->cd.iatmr1);
    /* also, the first N in the chain is different: */
    for (i = 0; i < iatmr1[0]; i++) {
      char *atmnam = get_coord_atmnam(&struc->cd, i);
      if (strcmp(atmnam, "N") == 0) {
        iatta[i] = iatmcls(0, "NH3", atclass, libs);
      }
      g_free(atmnam);
    }

    /* also, the O in the last residue are different: */
    for (i = iatmr1[seq->nres - 1]; i < struc->cd.natm; i++) {
      char *atmnam = get_coord_atmnam(&struc->cd, i);
      if (strcmp(atmnam, "OT") == 0 || strcmp(atmnam, "OT1") == 0
          || strcmp(atmnam, "OT2") == 0 || strcmp(atmnam, "OXT") == 0
          || strcmp(atmnam, "O") == 0) {
        iatta[i] = iatmcls(0, "OT1", atclass, libs);
      }
      g_free(atmnam);
    }
  }
}

static int *make_atom_type(const struct alignment *aln, int is,
                                 const struct mdt_library *mlib,
                                 const struct mdt_atom_class_list *atclass,
                                 int ifi, const struct libraries *libs)
{
  int *iatta;
  struct structure *struc = alignment_structure_get(aln, is);
  struct sequence *seq = alignment_sequence_get(aln, is);
  iatta = g_malloc(sizeof(int) * struc->cd.natm);
  atmcls_special(struc, seq, iatta, atclass, mlib, libs);
  return iatta;
}

/** Get/calculate the array of atom type bin indices */
static const int *property_iatta(const struct alignment *aln, int is,
                                 struct mdt_properties *prop,
                                 const struct mdt_library *mlib, int ifi,
                                 const struct libraries *libs)
{
  is--;
  if (!prop[is].iatta) {
    prop[is].iatta = make_atom_type(aln, is, mlib, mlib->atclass[0], ifi, libs);
  }
  return prop[is].iatta;
}

/** Get/calculate the array of hydrogen bond atom type bin indices */
static const int *property_hb_iatta(const struct alignment *aln, int is,
                                    struct mdt_properties *prop,
                                    const struct mdt_library *mlib, int ifi,
                                    const struct libraries *libs)
{
  is--;
  if (!prop[is].hb_iatta) {
    prop[is].hb_iatta = make_atom_type(aln, is, mlib, mlib->hbond, ifi, libs);
  }
  return prop[is].hb_iatta;
}

/** Get/calculate the hydrogen bond satisfaction index */
static float property_hbpot(const struct alignment *aln, int is,
                            struct mdt_properties *prop,
                            const struct mdt_library *mlib, int ifi,
                            const struct libraries *libs)
{
  struct structure *struc = alignment_structure_get(aln, is);
  const int *iatta = property_hb_iatta(aln, is, prop, mlib, ifi, libs);
  is--;
  if (!prop[is].hbpot) {
    prop[is].hbpot = g_malloc(sizeof(float));
    *(prop[is].hbpot) = hb_satisfaction(&struc->cd, iatta, mlib->hbond,
                                        mlib->hbond_cutoff);
  }
  return *(prop[is].hbpot);
}

/** Get/calculate the array of atom accessibility bin indices */
static const int *property_iatmacc(const struct alignment *aln, int is,
                                   struct mdt_properties *prop,
                                   const struct mdt_library *mlib, int ifi,
                                   const struct mdt_feature *feat)
{
  is--;
  if (!prop[is].iatmacc) {
    struct structure *struc = alignment_structure_get(aln, is);
    prop[is].iatmacc = g_malloc(sizeof(int) * struc->cd.natm);
    alliclsbin(struc->cd.natm, f_float1_pt(&struc->cd.atmacc), prop[is].iatmacc,
               mlib, ifi, feat->nbins - 1);
  }
  return prop[is].iatmacc;
}

/** Get/calculate the array of fractional atom accessibility bin indices */
static const int *property_ifatmacc(const struct alignment *aln, int is,
                                    struct mdt_properties *prop,
                                    const struct mdt_library *mlib, int ifi,
                                    const struct mdt_feature *feat,
                                    const struct libraries *libs, GError **err)
{
  is--;
  if (!prop[is].ifatmacc) {
    int i, *ifatmacc;
    struct sequence *seq = alignment_sequence_get(aln, is);
    struct structure *struc = alignment_structure_get(aln, is);

    ifatmacc = g_malloc(sizeof(int) * struc->cd.natm);
    for (i = 0; i < struc->cd.natm; i++) {
      int iattyp, ierr;
      float r, fatmacc;

      /* Get integer atom type */
      iattyp = coordinates_atom_type_get(&struc->cd, seq, i, libs, &ierr);
      if (ierr) {
        handle_modeller_error(err);
        g_free(ifatmacc);
        return NULL;
      }
      /* Get VDW atom radius */
      r = f_float2_get(&libs->vdwcnt, iattyp - 1, libs->tpl.submodel - 1);
      /* Calculate fractional atom accessibility from raw values */
      fatmacc = f_float1_get(&struc->cd.atmacc, i) / (4. * G_PI * r * r);
      /* Get the corresponding bin index */
      ifatmacc[i] = iclsbin(fatmacc, mlib, ifi, feat->nbins - 1);
    }
    prop[is].ifatmacc = ifatmacc;
  }
  return prop[is].ifatmacc;
}

static float dist1(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float xd, yd, zd;
  xd = x1 - x2;
  yd = y1 - y2;
  zd = z1 - z2;
  return sqrt(xd * xd + yd * yd + zd * zd);
}

/** Return the bin index for the distance between two specified atoms in the
    same protein. */
static int idist0(int ia1, int ia1p, const struct structure *struc,
                  const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 > 0 && ia1p > 0) {
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    ia1--;
    ia1p--;
    d = dist1(x[ia1], y[ia1], z[ia1], x[ia1p], y[ia1p], z[ia1p]);
    return iclsbin(d, mlib, ifi, nrang);
  } else {
    return nrang + 1;
  }
}

int my_mdt_index(int ifi, const struct alignment *aln, int is1, int ip1,
                 int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                 int ia1p, const struct mdt_library *mlib, int ip2,
                 int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                 const struct libraries *libs,
                 const struct energy_data *edat,
                 struct mdt_properties *prop, GError **err)
{
  int ret, ierr = 0;
  const int *binprop;
  struct structure *struc1, *struc2;
  struct sequence *seq1, *seq2;
  struct mdt_feature *feat = &mlib->base.features[ifi-1];
  struc1 = alignment_structure_get(aln, is1-1);
  struc2 = alignment_structure_get(aln, is2-1);
  seq1 = alignment_sequence_get(aln, is1-1);
  seq2 = alignment_sequence_get(aln, is2-1);
  switch(ifi) {
  case 66:
    return irestab(&aln->ialn, aln->naln, is1, f_int1_pt(&seq1->irestyp),
                   seq1->nres, ip1, mlib->deltai, mlib->deltai_ali,
                   feat->nbins, libs->igaptyp);
  case 67:
    return irestab(&aln->ialn, aln->naln, is2, f_int1_pt(&seq2->irestyp),
                   seq2->nres, ip1, mlib->deltai, mlib->deltai_ali,
                   feat->nbins, libs->igaptyp);
  case 77:
    return irestab(&aln->ialn, aln->naln, is1, f_int1_pt(&seq1->irestyp),
                   seq1->nres, ip1, mlib->deltaj, mlib->deltaj_ali,
                   feat->nbins, libs->igaptyp);
  case 78:
    return irestab(&aln->ialn, aln->naln, is2, f_int1_pt(&seq2->irestyp),
                   seq2->nres, ip1, mlib->deltaj, mlib->deltaj_ali,
                   feat->nbins, libs->igaptyp);
  case 79:
    return itable(property_iatta(aln, is1, prop, mlib, ifi, libs),
                  struc1->cd.natm, ia1, feat->nbins);
  case 80:
    return itable(property_iatmacc(aln, is1, prop, mlib, ifi, feat),
                  struc1->cd.natm, ia1, feat->nbins);
  case 81:
    binprop = property_ifatmacc(aln, is1, prop, mlib, ifi, feat, libs, err);
    if (!binprop) {
      return 0.;
    }
    return itable(binprop, struc1->cd.natm, ia1, feat->nbins);
  case 82: case 103:
    return idist0(ia1, ia1p, struc1, mlib, ifi, feat->nbins);
  case 83:
    return itable(property_iatta(aln, is1, prop, mlib, ifi, libs),
                  struc1->cd.natm, ia1p, feat->nbins);
  case 84:
    return numb_hda(ia1, property_hb_iatta(aln, is1, prop, mlib, ifi, libs),
                    &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    0, feat->nbins);
  case 85:
    return numb_hda(ia1, property_hb_iatta(aln, is1, prop, mlib, ifi, libs),
                    &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    1, feat->nbins);
  case 86:
    return iclsbin(property_hbpot(aln, is1, prop, mlib, ifi, libs), mlib, ifi,
                   feat->nbins - 1);
  case 87:
    return numb_hda(ia1, property_hb_iatta(aln, is1, prop, mlib, ifi, libs),
                    &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    2, feat->nbins);
  case 93: case 95: case 97: case 99:
    return itable(f_int1_pt(&struc1->iacc), seq1->nres, ir1, feat->nbins);
  case 94: case 96: case 98: case 100:
    return itable(f_int1_pt(&struc2->iacc), seq2->nres, ir2, feat->nbins);
  default:
    ret = mdt_index(ifi, aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1, ia1p,
                    &mlib->base, ip2, ibnd1, ibnd1p, is3, ir3, ir3p, libs, edat,
                    &ierr);
    if (ierr) {
      handle_modeller_error(err);
    }
    return ret;
  }
}
