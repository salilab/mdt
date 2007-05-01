/** \file mdt_property.c  Functions to precalculate protein properties used
 *                        to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <string.h>
#include <math.h>
#include "modeller.h"
#include "util.h"
#include "mdt_index.h"
#include "mdt_property.h"
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
    prop[i].iresol = 0;
    prop[i].radius_gyration = -1;
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
                           int ifi, const struct libraries *libs, GError **err)
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
const int *property_iatta(const struct alignment *aln, int is,
                          struct mdt_properties *prop,
                          const struct mdt_library *mlib, int ifi,
                          const struct libraries *libs, GError **err)
{
  if (!prop[is].iatta) {
    prop[is].iatta = make_atom_type(aln, is, mlib, mlib->atclass[0], ifi, libs,
                                    err);
  }
  return prop[is].iatta;
}

/** Get/calculate the array of hydrogen bond atom type bin indices */
const int *property_hb_iatta(const struct alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct libraries *libs, GError **err)
{
  if (!prop[is].hb_iatta) {
    prop[is].hb_iatta = make_atom_type(aln, is, mlib, mlib->hbond, ifi, libs,
                                       err);
  }
  return prop[is].hb_iatta;
}

/** Get/calculate the hydrogen bond satisfaction index */
gboolean property_hbpot(const struct alignment *aln, int is,
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
  if (!prop[is].hbpot) {
    prop[is].hbpot = g_malloc(sizeof(float));
    *(prop[is].hbpot) = hb_satisfaction(&struc->cd, iatta, mlib->hbond,
                                        mlib->hbond_cutoff);
  }
  *hbpot = *(prop[is].hbpot);
  return TRUE;
}

/** Get/calculate the resolution bin index */
int property_iresol(const struct alignment *aln, int is,
                    struct mdt_properties *prop,
                    const struct mdt_library *mlib, int ifi,
                    const struct mdt_libfeature *feat)
{
  if (prop[is].iresol == 0) {
    struct sequence *seq = alignment_sequence_get(aln, is);
    float resol;
    int iresol;

    /* artificially change the resolution of the NMR structures
       from the defined -1.00 to 0.45, to decrease the number of
       bins required to hold all defined resolutions while still
       separating NMR from X-ray structures: */
    resol = (seq->resol == -1.00 ? 0.45 : seq->resol);

    alliclsbin(1, &resol, &iresol, mlib, ifi, feat->nbins - 1);
    prop[is].iresol = iresol;
  }
  return prop[is].iresol;
}

/** Get center of mass */
static void get_mass_center(const float x[], const float y[], const float z[],
                            int natm, float *cx, float *cy, float *cz)
{
  int i;
  *cx = *cy = *cz = 0.;
  for (i = 0; i < natm; i++) {
    *cx += x[i];
    *cy += y[i];
    *cz += z[i];
  }
  if (natm > 0) {
    *cx /= natm;
    *cy /= natm;
    *cz /= natm;
  }
}

/** Get radius of gyration */
static float get_radius_gyration(const float x[], const float y[],
                                 const float z[], int natm, float cx,
                                 float cy, float cz)
{
  int i;
  float sum = 0.;
  for (i = 0; i < natm; i++) {
    float dx, dy, dz;
    dx = x[i] - cx;
    dy = y[i] - cy;
    dz = z[i] - cz;
    sum += dx * dx + dy * dy + dz * dz;
  }
  if (natm > 0) {
    sum /= natm;
  }
  return sqrt(sum);
}

/** Get/calculate the radius of gyration bin index */
int property_radius_gyration(const struct alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct mdt_libfeature *feat)
{
  if (prop[is].radius_gyration == -1) {
    struct structure *struc = alignment_structure_get(aln, is);
    float *x, *y, *z, radius_gyration, cx, cy, cz;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    /* get center of mass */
    get_mass_center(x, y, z, struc->cd.natm, &cx, &cy, &cz);
    /* get radius of gyration */
    radius_gyration = get_radius_gyration(x, y, z, struc->cd.natm, cx, cy, cz);
    alliclsbin(1, &radius_gyration, &prop[is].radius_gyration, mlib, ifi,
               feat->nbins - 1);
  }
  return prop[is].radius_gyration;
}

/** Get/calculate the array of atom accessibility bin indices */
const int *property_iatmacc(const struct alignment *aln, int is,
                            struct mdt_properties *prop,
                            const struct mdt_library *mlib, int ifi,
                            const struct mdt_libfeature *feat)
{
  if (!prop[is].iatmacc) {
    struct structure *struc = alignment_structure_get(aln, is);
    prop[is].iatmacc = g_malloc(sizeof(int) * struc->cd.natm);
    alliclsbin(struc->cd.natm, f_float1_pt(&struc->cd.atmacc),
               prop[is].iatmacc, mlib, ifi, feat->nbins - 1);
  }
  return prop[is].iatmacc;
}

/** Get/calculate the array of fractional atom accessibility bin indices */
const int *property_ifatmacc(const struct alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib, int ifi,
                             const struct mdt_libfeature *feat,
                             const struct libraries *libs, GError **err)
{
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
const struct mdt_bond *property_one_bond(const struct alignment *aln,
                                         int is, struct mdt_properties *prop,
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
                                                 const struct mdt_library
                                                 *mlib,
                                                 const struct libraries *libs)
{
  if (!prop[is].triplets) {
    struct sequence *seq = alignment_sequence_get(aln, is);
    struct structure *struc = alignment_structure_get(aln, is);
    prop[is].triplets = trpclass(struc, seq, mlib->trpclass, libs);
  }
  return prop[is].triplets;
}

/** Require that the triplets have at least min_natom atoms each */
gboolean triplet_require_natom(const struct mdt_library *mlib, int min_natom,
                               int ifeat, GError **err)
{
  int have_natom = mlib->trpclass->natom;
  if (have_natom < min_natom) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "This feature (%d) works only for tuplets of %d or more "
                "atoms; you have only %d", ifeat, min_natom, have_natom);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Get a single atom triplet from a structure */
const struct mdt_triplet *property_one_triplet(const struct alignment *aln,
                                               int is,
                                               struct mdt_properties *prop,
                                               const struct mdt_library *mlib,
                                               int ibnd1, int ia1,
                                               const struct libraries *libs)
{
  const struct mdt_triplet_list *trp;
  trp = property_triplets(aln, is, prop, mlib, libs);

  return &trp[ia1].triplets[ibnd1];
}
