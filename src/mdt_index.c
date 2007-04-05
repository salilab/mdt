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
#include "mdt_stereo.h"
#include "mdt_triplets.h"

/** Make a new mdt_properties structure */
struct mdt_properties *mdt_properties_new(const struct alignment *aln)
{
  struct mdt_properties *prop;
  int i, j;
  prop = g_malloc(sizeof(struct mdt_properties) * aln->naln);
  for (i = 0; i < aln->naln; i++) {
    for (j = 0; j < N_MDT_BOND_TYPES; j++) {
      prop[i].bonds[j] = NULL;
    }
    prop[i].triplets = NULL;
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
  int i, j;
  for (i = 0; i < aln->naln; i++) {
    struct structure *struc = alignment_structure_get(aln, i);
    for (j = 0; j < N_MDT_BOND_TYPES; j++) {
      if (prop[i].bonds[j]) {
        g_free(prop[i].bonds[j]->bonds);
      }
      g_free(prop[i].bonds[j]);
    }
    if (prop[i].triplets) {
      for (j = 0; j < struc->cd.natm; j++) {
        g_free(prop[i].triplets[j].triplets);
      }
    }
    g_free(prop[i].triplets);
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

static void atmclass_disulfide(const int iss[], int nss,
                               const struct structure *struc,
                               const struct sequence *seq,
                               const struct mdt_atom_class_list *atclass,
                               int iatta[], const struct libraries *libs)
{
  int cycint, i, *iatmr1;

  cycint = residue_type_from_name("CSS", libs);

  iatmr1 = f_int1_pt(&struc->cd.iatmr1);
  for (i = 0; i < nss; i++) {
    int ir1;
    for (ir1 = 0; ir1 < 2; i++) {
      int iatm, istart, iend, ir = iss[i * 2 + ir1];
      istart = iatmr1[ir - 1] - 1;
      if (ir < seq->nres) {
        iend = iatmr1[ir];
      } else {
        iend = struc->cd.natm;
      }
      for (iatm = istart; iatm < iend; iatm++) {
        char *atmnam = get_coord_atmnam(&struc->cd, iatm);
        iatta[iatm] = iatmcls(cycint, atmnam, atclass, libs);
        g_free(atmnam);
      }
    }
  }
}

static gboolean atmcls_special(struct structure *struc,
                               const struct sequence *seq, int iatta[],
                               const struct mdt_atom_class_list *atclass,
                               const struct mdt_library *mlib,
                               const struct libraries *libs, GError **err)
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

  if (mlib->special_atoms) {
    int *iss, nss, ierr, *iatmr1 = f_int1_pt(&struc->cd.iatmr1);
    /* Take care of the atoms in the disulfide bonded Cys residues: */
    mod_find_ss(&iss, &nss, struc, seq, &ierr);
    if (ierr != 0) {
      handle_modeller_error(err);
      return FALSE;
    }
    if (nss > 0) {
      atmclass_disulfide(iss, nss, struc, seq, atclass, iatta, libs);
      g_free(iss);
    }

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
  return TRUE;
}

static int *make_atom_type(const struct alignment *aln, int is,
                                 const struct mdt_library *mlib,
                                 const struct mdt_atom_class_list *atclass,
                                 int ifi, const struct libraries *libs,
                                 GError **err)
{
  int *iatta;
  struct structure *struc = alignment_structure_get(aln, is);
  struct sequence *seq = alignment_sequence_get(aln, is);
  iatta = g_malloc(sizeof(int) * struc->cd.natm);
  if (!atmcls_special(struc, seq, iatta, atclass, mlib, libs, err)) {
    g_free(iatta);
    return NULL;
  } else {
    return iatta;
  }
}

/** Get/calculate the array of atom type bin indices */
static const int *property_iatta(const struct alignment *aln, int is,
                                 struct mdt_properties *prop,
                                 const struct mdt_library *mlib, int ifi,
                                 const struct libraries *libs, GError **err)
{
  is--;
  if (!prop[is].iatta) {
    prop[is].iatta = make_atom_type(aln, is, mlib, mlib->atclass[0], ifi, libs,
                                    err);
  }
  return prop[is].iatta;
}

/** Get/calculate the array of hydrogen bond atom type bin indices */
static const int *property_hb_iatta(const struct alignment *aln, int is,
                                    struct mdt_properties *prop,
                                    const struct mdt_library *mlib, int ifi,
                                    const struct libraries *libs, GError **err)
{
  is--;
  if (!prop[is].hb_iatta) {
    prop[is].hb_iatta = make_atom_type(aln, is, mlib, mlib->hbond, ifi, libs,
                                       err);
  }
  return prop[is].hb_iatta;
}

/** Get/calculate the hydrogen bond satisfaction index */
static gboolean property_hbpot(const struct alignment *aln, int is,
                               struct mdt_properties *prop,
                               const struct mdt_library *mlib, int ifi,
                               const struct libraries *libs, float *hbpot,
                               GError **err)
{
  struct structure *struc = alignment_structure_get(aln, is);
  const int *iatta = property_hb_iatta(aln, is, prop, mlib, ifi, libs, err);
  if (!iatta) {
    return FALSE;
  }
  is--;
  if (!prop[is].hbpot) {
    prop[is].hbpot = g_malloc(sizeof(float));
    *(prop[is].hbpot) = hb_satisfaction(&struc->cd, iatta, mlib->hbond,
                                        mlib->hbond_cutoff);
  }
  *hbpot = *(prop[is].hbpot);
  return TRUE;
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

/** Get/calculate the list of all bonds for a structure. */
const struct mdt_bond_list *property_bonds(const struct alignment *aln, int is,
                                           struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int bondtype,
                                           const struct libraries *libs)
{
  is--;
  if (!prop[is].bonds[bondtype]) {
    struct sequence *seq = alignment_sequence_get(aln, is);
    struct structure *struc = alignment_structure_get(aln, is);
    prop[is].bonds[bondtype] = get_stereo(struc, seq,
                                          mlib->atclass[bondtype + 1],
                                          bondtype, libs);
  }
  return prop[is].bonds[bondtype];
}

/** Get a single bond from a structure */
static const struct mdt_bond *property_one_bond(const struct alignment *aln,
                                                int is,
                                                struct mdt_properties *prop,
                                                const struct mdt_library *mlib,
                                                int bondtype, int ibnd1,
                                                const struct libraries *libs)
{
  return &property_bonds(aln, is, prop, mlib, bondtype, libs)->bonds[ibnd1];
}


/** Get/calculate the list of all triplets for a structure. */
const struct mdt_triplet_list *property_triplets(const struct alignment *aln,
                                                 int is,
                                                 struct mdt_properties *prop,
                                                 const struct mdt_library *mlib,
                                                 const struct libraries *libs)
{
  is--;
  if (!prop[is].triplets) {
    struct sequence *seq = alignment_sequence_get(aln, is);
    struct structure *struc = alignment_structure_get(aln, is);
    prop[is].triplets = trpclass(struc, seq, mlib->trpclass, libs);
  }
  return prop[is].triplets;
}

/** Get a single atom triplet from a structure */
static const struct mdt_triplet
    *property_one_triplet(const struct alignment *aln, int is,
                          struct mdt_properties *prop,
                          const struct mdt_library *mlib, int ibnd1, int ia1,
                          const struct libraries *libs)
{
  const struct mdt_triplet_list *trp;
  trp = property_triplets(aln, is, prop, mlib, libs);

  return &trp[ia1-1].triplets[ibnd1-1];
}

static float dist1(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float xd, yd, zd;
  xd = x1 - x2;
  yd = y1 - y2;
  zd = z1 - z2;
  return sqrt(xd * xd + yd * yd + zd * zd);
}

static float angle1(float x1, float y1, float z1, float x2, float y2, float z2,
                    float x3, float y3, float z3)
{
  static const float tiny = 1.0e-15;
  float d1, d2, v1x, v1y, v1z, v2x, v2y, v2z, scalprod, sizeprod, div;
  v1x = x1-x2;
  v1y = y1-y2;
  v1z = z1-z2;
  v2x = x3-x2;
  v2y = y3-y2;
  v2z = z3-z2;
  d1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
  d2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
  scalprod = v1x*v2x+v1y*v2y+v1z*v2z;
  sizeprod = d1 * d2;
  div = (sizeprod > tiny ? scalprod / sizeprod : 0.0);
  div = CLAMP(div, -1.0, 1.0);
  return acos(div) * 180.0 / G_PI;
}

static float dihedral1(float x1, float y1, float z1, float x2, float y2,
                       float z2, float x3, float y3, float z3, float x4,
                       float y4, float z4)
{
  double rt[4][3], l1[3], l2[3], l3[3], xt1[3], xt2[3], leng1, leng2, dot1, ang,
         sign, norm;

  rt[0][0]=x1;
  rt[0][1]=y1;
  rt[0][2]=z1;
  rt[1][0]=x2;
  rt[1][1]=y2;
  rt[1][2]=z2;
  rt[2][0]=x3;
  rt[2][1]=y3;
  rt[2][2]=z3;
  rt[3][0]=x4;
  rt[3][1]=y4;
  rt[3][2]=z4;

  l1[0]=rt[1][0]-rt[0][0];
  l1[1]=rt[1][1]-rt[0][1];
  l1[2]=rt[1][2]-rt[0][2];
  l2[0]=rt[2][0]-rt[1][0];
  l2[1]=rt[2][1]-rt[1][1];
  l2[2]=rt[2][2]-rt[1][2];
  l3[0]=rt[3][0]-rt[2][0];
  l3[1]=rt[3][1]-rt[2][1];
  l3[2]=rt[3][2]-rt[2][2];

  xt1[0]=l2[1]*l1[2]-l2[2]*l1[1];
  xt1[1]=l2[2]*l1[0]-l2[0]*l1[2];
  xt1[2]=l2[0]*l1[1]-l2[1]*l1[0];
  xt2[0]=l3[1]*l2[2]-l3[2]*l2[1];
  xt2[1]=l3[2]*l2[0]-l3[0]*l2[2];
  xt2[2]=l3[0]*l2[1]-l3[1]*l2[0];
  leng1 = xt1[0] * xt1[0] + xt1[1] * xt1[1] + xt1[2] * xt1[2];
  leng2 = xt2[0] * xt2[0] + xt2[1] * xt2[1] + xt2[2] * xt2[2];
  dot1 = xt1[0]*xt2[0]+xt1[1]*xt2[1]+xt1[2]*xt2[2];
  norm = dot1/sqrt(leng1*leng2);
  norm = CLAMP(norm, -1.0, 1.0);
  ang = acos(norm);
  sign=xt1[0]*l3[0]+xt1[1]*l3[1]+xt1[2]*l3[2];
  if (sign < 0.0) {
    ang = -ang;
  }
  return -ang * 180.0 / G_PI;
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

/** Return the bin index for the angle between three specified atoms in the
    same protein. */
static int iangle0(int ia1, int ia2, int ia3, const struct structure *struc,
                   const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0) {
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    d = angle1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
               y[ia3], z[ia3]);
    return iclsbin(d, mlib, ifi, nrang);
  } else {
    return nrang + 1;
  }
}

/** Return the bin index for the dihedral angle between four specified atoms
    in the same protein. */
static int idihedral0(int ia1, int ia2, int ia3, int ia4,
                      const struct structure *struc,
                      const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0 && ia4 >= 0) {
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    d = dihedral1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
                  y[ia3], z[ia3], x[ia4], y[ia4], z[ia4]);
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
  float fprop;
  struct structure *struc1, *struc2;
  struct sequence *seq1, *seq2;
  const struct mdt_bond *bond;
  const struct mdt_triplet *trp, *trp2;
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
    binprop = property_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0.;
    }
    return itable(binprop, struc1->cd.natm, ia1, feat->nbins);
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
    binprop = property_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0.;
    }
    return itable(binprop, struc1->cd.natm, ia1p, feat->nbins);
  case 84:
    binprop = property_hb_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0.;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    0, feat->nbins);
  case 85:
    binprop = property_hb_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0.;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    1, feat->nbins);
  case 86:
    if (!property_hbpot(aln, is1, prop, mlib, ifi, libs, &fprop, err)) {
      return 0.;
    }
    return iclsbin(fprop, mlib, ifi, feat->nbins - 1);
  case 87:
    binprop = property_hb_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0.;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    2, feat->nbins);
  case 93: case 95: case 97: case 99:
    return itable(f_int1_pt(&struc1->iacc), seq1->nres, ir1, feat->nbins);
  case 94: case 96: case 98: case 100:
    return itable(f_int1_pt(&struc2->iacc), seq2->nres, ir2, feat->nbins);
  case 101:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return CLAMP(trp->trpclass, 1, feat->nbins);
  case 102:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return CLAMP(trp->trpclass, 1, feat->nbins);
  case 104:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return iangle0(ia1-1, ia1p-1, trp->iata[0], struc1, mlib, ifi, feat->nbins);
  case 105:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return iangle0(trp->iata[0], ia1-1, ia1p-1, struc1, mlib, ifi, feat->nbins);
  case 106:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    trp2 = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(trp->iata[0], ia1-1, ia1p-1, trp2->iata[0], struc1,
                      mlib, ifi, feat->nbins);
  case 107:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    trp2 = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(trp->iata[1], trp2->iata[0], ia1-1, ia1p-1, struc1,
                      mlib, ifi, feat->nbins);
  case 108:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    trp2 = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(ia1-1, ia1p-1, trp->iata[0], trp2->iata[1], struc1,
                      mlib, ifi, feat->nbins);
  case 109:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return CLAMP(bond->bndgrp, 1, feat->nbins);
  case 110:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return idist0(bond->iata[0] + 1, bond->iata[1] + 1, struc1, mlib, ifi,
                  feat->nbins);
  case 111:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return CLAMP(bond->bndgrp, 1, feat->nbins);
  case 112:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return iangle0(bond->iata[0], bond->iata[1], bond->iata[2], struc1, mlib,
                   ifi, feat->nbins);
  case 113:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return CLAMP(bond->bndgrp, 1, feat->nbins);
  case 114:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return idihedral0(bond->iata[0], bond->iata[1], bond->iata[2],
                      bond->iata[3], struc1, mlib, ifi, feat->nbins);
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
