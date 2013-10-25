/** \file mdt_property.c  Functions to precalculate protein properties used
 *                        to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <string.h>
#include <math.h>
#include "modeller.h"
#include "util.h"
#include "mdt_index.h"
#include "mdt_disulfides.h"
#include "mdt_property.h"
#include "mdt_hydrogen_bonds.h"
#include "mdt_stereo.h"
#include "mdt_tuples.h"
#include "mdt_atom_classes.h"

/** Precalculated per-sequence properties for calculating MDT indices */
struct mdt_properties {
  /** Lists of bonds */
  struct mdt_bond_list *bonds[N_MDT_BOND_TYPES];
  /** Hash of excluded atom pairs */
  GHashTable *exclusions;
  /** Callbacks for writing properties to HDF5 files */
  GHashTable *write_lib_funcs;
  /** Lists of atom tuples for each atom */
  struct mdt_tuple_list *tuples;
  /** Bin indices for hydrogen bond atom type */
  int *hb_iatta;
  /** Hydrogen bond satisfaction index */
  float *hbpot;
  /** Radius of gyration, or -1 if not yet calculated */
  float radius_gyration;
  /** Bin indices for atom type */
  int *iatta;
  /** Fractional atom accessibility */
  float *fatmacc;
  /** Average sidechain Biso */
  float *sidechain_biso;
  /** Atom indices for distance_atoms */
  int *dstind1, *dstind2;
  /** Atom types for residue bond separation */
  int *resbond_attyp;
  /** Residue indices of S-S bonds */
  struct mdt_disulfide_list *disulfides;
  /** Number of user-defined properties */
  int n_user_properties;
  /** User-defined properties */
  float **user_properties;
};

/** Make a new mdt_properties structure */
struct mdt_properties *mdt_properties_new(const struct mod_alignment *aln,
                                          const struct mdt_library *mlib)
{
  struct mdt_properties *prop;
  int i, j;
  prop = g_malloc(sizeof(struct mdt_properties) * aln->nseq);
  for (i = 0; i < aln->nseq; i++) {
    for (j = 0; j < N_MDT_BOND_TYPES; j++) {
      prop[i].bonds[j] = NULL;
    }
    prop[i].exclusions = NULL;
    prop[i].write_lib_funcs = g_hash_table_new(NULL, NULL);
    prop[i].tuples = NULL;
    prop[i].hb_iatta = NULL;
    prop[i].hbpot = NULL;
    prop[i].radius_gyration = -1;
    prop[i].iatta = NULL;
    prop[i].fatmacc = NULL;
    prop[i].sidechain_biso = NULL;
    prop[i].dstind1 = NULL;
    prop[i].dstind2 = NULL;
    prop[i].resbond_attyp = NULL;
    prop[i].disulfides = NULL;
    prop[i].n_user_properties = mlib->user_properties->len;
    prop[i].user_properties = g_malloc0(sizeof(float *)
                                        * prop[i].n_user_properties);
  }
  return prop;
}

/** Free an mdt_properties structure */
void mdt_properties_free(struct mdt_properties *prop,
                         const struct mod_alignment *aln)
{
  int i, j;
  for (i = 0; i < aln->nseq; i++) {
    struct mod_structure *struc = mod_alignment_structure_get(aln, i);
    for (j = 0; j < N_MDT_BOND_TYPES; j++) {
      if (prop[i].bonds[j]) {
        g_free(prop[i].bonds[j]->bonds);
      }
      g_free(prop[i].bonds[j]);
    }
    if (prop[i].exclusions) {
      g_hash_table_destroy(prop[i].exclusions);
    }
    if (prop[i].write_lib_funcs) {
      g_hash_table_destroy(prop[i].write_lib_funcs);
    }
    if (prop[i].tuples) {
      for (j = 0; j < struc->cd.natm; j++) {
        g_free(prop[i].tuples[j].tuples);
      }
    }
    if (prop[i].disulfides) {
      g_free(prop[i].disulfides->ss);
    }
    g_free(prop[i].tuples);
    g_free(prop[i].hb_iatta);
    g_free(prop[i].hbpot);
    g_free(prop[i].iatta);
    g_free(prop[i].fatmacc);
    g_free(prop[i].sidechain_biso);
    g_free(prop[i].dstind1);
    g_free(prop[i].dstind2);
    g_free(prop[i].resbond_attyp);
    g_free(prop[i].disulfides);
    for (j = 0; j < prop[i].n_user_properties; ++j) {
      g_free(prop[i].user_properties[j]);
    }
    g_free(prop[i].user_properties);
  }
  g_free(prop);
}

static int iatmcls(int irestyp, const char *atmnam,
                   const struct mdt_atom_class_list *atclass,
                   const struct mod_libraries *libs)
{
  gboolean allres, allatm;
  char *resnam;
  int iclass;

  allres = (irestyp == 0);
  allatm = (strcmp(atmnam, "*") == 0);

  if (allres) {
    resnam = g_strdup("*");
  } else {
    resnam = mod_residue_name_from_type(irestyp, libs);
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

  if (!mod_residue_is_hetatm(irestyp, libs) && !mod_atom_is_hydrogen(atmnam)) {
    mod_logwarning("iatmcls", "Model atom not classified: %s:%s", resnam,
                   atmnam);
  }
  g_free(resnam);
  return 0;
}

static void atmclass_disulfide(const int iss[], int nss,
                               const struct mod_structure *struc,
                               const struct mod_sequence *seq,
                               const struct mdt_atom_class_list *atclass,
                               int iatta[], const struct mod_libraries *libs)
{
  int cycint, i, *iatmr1;

  cycint = mod_residue_type_from_name("CSS", libs);

  iatmr1 = mod_int1_pt(&struc->cd.iatmr1);
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
        char *atmnam = mod_coordinates_atmnam_get(&struc->cd, iatm);
        iatta[iatm] = iatmcls(cycint, atmnam, atclass, libs);
        g_free(atmnam);
      }
    }
  }
}

static gboolean atmcls_special(struct mod_structure *struc,
                               const struct mod_sequence *seq, int iatta[],
                               const struct mdt_atom_class_list *atclass,
                               const struct mdt_library *mlib,
                               const struct mod_libraries *libs, GError **err)
{
  int i, *irestyp, *iresatm;
  iresatm = mod_int1_pt(&struc->cd.iresatm);
  irestyp = mod_int1_pt(&seq->irestyp);
  for (i = 0; i < struc->cd.natm; i++) {
    int irest = irestyp[iresatm[i] - 1];
    char *atmnam = mod_coordinates_atmnam_get(&struc->cd, i);
    iatta[i] = iatmcls(irest, atmnam, atclass, libs);
    g_free(atmnam);
  }

  if (mlib->special_atoms) {
    int *iss, nss, ierr, *iatmr1 = mod_int1_pt(&struc->cd.iatmr1);
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
      char *atmnam = mod_coordinates_atmnam_get(&struc->cd, i);
      if (strcmp(atmnam, "N") == 0) {
        iatta[i] = iatmcls(0, "NH3", atclass, libs);
      }
      g_free(atmnam);
    }

    /* also, the O in the last residue are different: */
    for (i = iatmr1[seq->nres - 1]; i < struc->cd.natm; i++) {
      char *atmnam = mod_coordinates_atmnam_get(&struc->cd, i);
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

static int *make_atom_type(const struct mod_alignment *aln, int is,
                           const struct mdt_library *mlib,
                           const struct mdt_atom_class_list *atclass,
                           const struct mod_libraries *libs, GError **err)
{
  int *iatta;
  struct mod_structure *struc = mod_alignment_structure_get(aln, is);
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
  iatta = g_malloc(sizeof(int) * struc->cd.natm);
  if (!atmcls_special(struc, seq, iatta, atclass, mlib, libs, err)) {
    g_free(iatta);
    return NULL;
  } else {
    return iatta;
  }
}

static void add_write_lib_callback(struct mdt_properties *prop, int is,
                                   mdt_cb_write_lib writelibfunc)
{
  g_hash_table_insert(prop[is].write_lib_funcs, writelibfunc,
                      GINT_TO_POINTER(1));
}

static void copy_callbacks(gpointer key, gpointer value, gpointer user_data)
{
  struct mdt *mdt = (struct mdt *)user_data;
  mdt_cb_write_lib writelibfunc = (mdt_cb_write_lib)key;
  mdt_set_write_lib_callback(mdt, writelibfunc);
}

void mdt_property_get_write_callbacks(const struct mdt_properties *prop,
                                      const struct mod_alignment *aln,
                                      struct mdt *mdt)
{
  int is;
  for (is = 0; is < aln->nseq; ++is) {
    g_hash_table_foreach(prop[is].write_lib_funcs, copy_callbacks, mdt);
  }
}

/** Get/calculate the array of atom type bin indices */
const int *property_iatta(const struct mod_alignment *aln, int is,
                          struct mdt_properties *prop,
                          const struct mdt_library *mlib,
                          const struct mod_libraries *libs, GError **err)
{
  if (!prop[is].iatta) {
    prop[is].iatta = make_atom_type(aln, is, mlib, mlib->atclass[0], libs, err);
    add_write_lib_callback(prop, is, mdt_atom_class_write);
  }
  return prop[is].iatta;
}

/** Get/calculate the array of hydrogen bond atom type bin indices */
const int *property_hb_iatta(const struct mod_alignment *aln, int is,
                             struct mdt_properties *prop,
                             const struct mdt_library *mlib,
                             const struct mod_libraries *libs, GError **err)
{
  if (!prop[is].hb_iatta) {
    prop[is].hb_iatta = make_atom_type(aln, is, mlib, mlib->hbond, libs, err);
    add_write_lib_callback(prop, is, mdt_hbond_write);
  }
  return prop[is].hb_iatta;
}

/** Get/calculate the hydrogen bond satisfaction index */
gboolean property_hbpot(const struct mod_alignment *aln, int is,
                        struct mdt_properties *prop,
                        const struct mdt_library *mlib,
                        const struct mod_libraries *libs, float *hbpot,
                        GError **err)
{
  struct mod_structure *struc = mod_alignment_structure_get(aln, is);
  const int *iatta = property_hb_iatta(aln, is, prop, mlib, libs, err);
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

/** Get/calculate the radius of gyration */
float property_radius_gyration(const struct mod_alignment *aln, int is,
                               struct mdt_properties *prop)
{
  if (prop[is].radius_gyration == -1) {
    struct mod_structure *struc = mod_alignment_structure_get(aln, is);
    float *x, *y, *z, cx, cy, cz;
    x = mod_float1_pt(&struc->cd.x);
    y = mod_float1_pt(&struc->cd.y);
    z = mod_float1_pt(&struc->cd.z);
    /* get center of mass */
    get_mass_center(x, y, z, struc->cd.natm, &cx, &cy, &cz);
    /* get radius of gyration */
    prop[is].radius_gyration = get_radius_gyration(x, y, z, struc->cd.natm,
                                                   cx, cy, cz);
  }
  return prop[is].radius_gyration;
}

/** Get/calculate the array of fractional atom accessibilities.
    \return NULL on failure. */
const float *property_fatmacc(const struct mod_alignment *aln, int is,
                              struct mdt_properties *prop,
                              const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *struc = mod_alignment_structure_get(aln, is);
  if (!prop[is].fatmacc) {
    float *fatmacc;
    int i;
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);

    fatmacc = g_malloc(sizeof(float) * struc->cd.natm);
    for (i = 0; i < struc->cd.natm; i++) {
      int iattyp, ierr;
      float r;

      /* Get integer atom type */
      iattyp = mod_coordinates_atom_type_get(&struc->cd, seq, i, libs, &ierr);
      if (ierr) {
        handle_modeller_error(err);
        g_free(fatmacc);
        return NULL;
      }
      /* Get VDW atom radius */
      r = mod_float2_get(&libs->vdwcnt, iattyp - 1, libs->tpl.submodel - 1);
      /* Calculate fractional atom accessibility from raw values */
      fatmacc[i] = mod_float1_get(&struc->cd.atmacc, i) / (4. * G_PI * r * r);
    }
    prop[is].fatmacc = fatmacc;
  }
  return prop[is].fatmacc;
}

/** Get/calculate the array of average sidechain Biso. */
const float *property_sidechain_biso(const struct mod_alignment *aln, int is,
                                     struct mdt_properties *prop)
{
  struct mod_structure *struc = mod_alignment_structure_get(aln, is);
  if (!prop[is].sidechain_biso) {
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
    float *biso;
    int i, natm, total_ok = 0;
    float total_biso = 0.;

    biso = g_malloc(sizeof(float) * seq->nres);
    for (i = 0; i < seq->nres; ++i) {
      mod_residue_sidechain_biso(&struc->cd, seq, i, &biso[i], &natm);
      if (natm > 0) {
        total_biso += biso[i];
        total_ok++;
      }
    }
    /* Heuristic: make sure the units are in the order of tens */
    if (total_ok > 0) {
      total_biso /= total_ok;
      if (total_biso < 2.) {
        for (i = 0; i < seq->nres; ++i) {
          biso[i] *= 4.0 * G_PI * G_PI;
        }
      }
    }
    prop[is].sidechain_biso = biso;
  }
  return prop[is].sidechain_biso;
}

/** Get/calculate the array of per-residue distance atom indices */
void property_distance_atom_indices(const struct mod_alignment *aln, int is,
                                    struct mdt_properties *prop,
                                    const struct mdt_library *mlib,
                                    const int **dstind1, const int **dstind2)
{
  if (!prop[is].dstind1) {
    struct mod_structure *struc = mod_alignment_structure_get(aln, is);
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
    int i, *d1, *d2;
    d1 = g_malloc(sizeof(int) * seq->nres);
    d2 = g_malloc(sizeof(int) * seq->nres);
    for (i = 0; i < seq->nres; ++i) {
      d1[i] = mod_residue_find_atom(&struc->cd, seq, i,
                                    mlib->distance_atoms[0]) - 1;
      d2[i] = mod_residue_find_atom(&struc->cd, seq, i,
                                    mlib->distance_atoms[1]) - 1;
    }
    prop[is].dstind1 = d1;
    prop[is].dstind2 = d2;
  }
  *dstind1 = prop[is].dstind1;
  *dstind2 = prop[is].dstind2;
}

/** Get/calculate the list of all bonds for a structure. */
const struct mdt_bond_list *property_bonds(const struct mod_alignment *aln,
                                           int is,
                                           struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int bondtype,
                                           const struct mod_libraries *libs)
{
  if (!prop[is].bonds[bondtype]) {
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
    struct mod_structure *struc = mod_alignment_structure_get(aln, is);
    prop[is].bonds[bondtype] = get_stereo(struc, seq,
                                          mlib->atclass[bondtype + 1],
                                          bondtype, libs);
    switch(bondtype) {
      case MDT_BOND_TYPE_BOND:
        add_write_lib_callback(prop, is, mdt_bond_class_write);
        break;
      case MDT_BOND_TYPE_ANGLE:
        add_write_lib_callback(prop, is, mdt_angle_class_write);
        break;
      case MDT_BOND_TYPE_DIHEDRAL:
        add_write_lib_callback(prop, is, mdt_dihedral_class_write);
        break;
    }
  }
  return prop[is].bonds[bondtype];
}

/** Get a single bond from a structure */
const struct mdt_bond *property_one_bond(const struct mod_alignment *aln,
                                         int is, struct mdt_properties *prop,
                                         const struct mdt_library *mlib,
                                         int bondtype, int ibnd1,
                                         const struct mod_libraries *libs)
{
  return &property_bonds(aln, is, prop, mlib, bondtype, libs)->bonds[ibnd1];
}

static void add_exclusions(GHashTable *h, const struct mod_alignment *aln,
                           int is, struct mdt_properties *prop,
                           const struct mdt_library *mlib,
                           int bondtype, int atind0, int atind1,
                           const struct mod_libraries *libs)
{
  int i;
  const struct mdt_bond_list *b = property_bonds(aln, is, prop, mlib,
                                                 bondtype, libs);
  for (i = 0; i < b->nbonds; ++i) {
    const struct mdt_bond *bond = &b->bonds[i];
    g_hash_table_insert(h, MAKE_HASH_KEY(bond->iata[atind0],
                                         bond->iata[atind1]),
                        GINT_TO_POINTER(1));
  }
}

/* Make sure that the hash key is large enough to store all possible
   atom indices */
static gboolean check_exclusion_hash_key_size(const struct mod_structure *s1,
                                              GError **err)
{
#if MDT_SIZEOF_POINTER == 8
  guint64 max_atoms = G_MAXUINT32;
#elif MDT_SIZEOF_POINTER == 4
  guint64 max_atoms = G_MAXUINT64;
#endif
  if ((guint64)s1->cd.natm > max_atoms) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "Too many atoms in protein (%d): cannot exclude atom pairs",
                s1->cd.natm);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Get/calculate all excluded atom pairs for a structure.
    NULL is returned if no pairs are excluded. */
gboolean property_exclusions(const struct mod_alignment *aln,
                             int is, struct mdt_properties *prop,
                             const struct mdt_library *mlib,
                             gboolean exclude_bonds,
                             gboolean exclude_angles,
                             gboolean exclude_dihedrals,
                             const struct mod_libraries *libs,
                             GHashTable **exclusions, GError **err)
{
  if (exclude_bonds || exclude_angles || exclude_dihedrals) {
    if (!prop[is].exclusions) {
      GHashTable *h;
      struct mod_structure *struc = mod_alignment_structure_get(aln, is);
      if (!check_exclusion_hash_key_size(struc, err)) {
        return FALSE;
      }

      h = g_hash_table_new_full(g_direct_hash, g_direct_equal, NULL, NULL);
      if (exclude_bonds) {
        add_exclusions(h, aln, is, prop, mlib, MDT_BOND_TYPE_BOND, 0, 1, libs);
      }
      if (exclude_angles) {
        add_exclusions(h, aln, is, prop, mlib, MDT_BOND_TYPE_ANGLE, 0, 2, libs);
      }
      if (exclude_dihedrals) {
        add_exclusions(h, aln, is, prop, mlib, MDT_BOND_TYPE_DIHEDRAL, 0, 3,
                       libs);
      }
      prop[is].exclusions = h;
    }
    *exclusions = prop[is].exclusions;
  } else {
    *exclusions = NULL;
  }
  return TRUE;
}

/** Get/calculate the list of all tuples for a structure. */
const struct mdt_tuple_list *property_tuples(const struct mod_alignment *aln,
                                             int is,
                                             struct mdt_properties *prop,
                                             const struct mdt_library *mlib,
                                             const struct mod_libraries *libs)
{
  if (!prop[is].tuples) {
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
    struct mod_structure *struc = mod_alignment_structure_get(aln, is);
    prop[is].tuples = tupclass(struc, seq, mlib->tupclass, libs);
    add_write_lib_callback(prop, is, mdt_tuple_write);
  }
  return prop[is].tuples;
}

/** Get/calculate the array of atom types for residue bond separation */
const int *property_resbond_attyp(const struct mod_alignment *aln, int is,
                                  struct mdt_properties *prop,
                                  const struct mdt_library *mlib,
                                  const struct mod_libraries *libs)
{
  if (!prop[is].resbond_attyp) {
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
    struct mod_structure *struc = mod_alignment_structure_get(aln, is);

    /* Populate the atom types (once per sequence) */
    prop[is].resbond_attyp = mdt_residue_bonds_assign_atom_types(struc, seq,
                                               &mlib->residue_bond_list, libs);
    add_write_lib_callback(prop, is, mdt_bond_class_write);
  }
  return prop[is].resbond_attyp;
}

/** Get/calculate the list of disulfide bridges */
const struct mdt_disulfide_list *property_disulfides(
                               const struct mod_alignment *aln, int is,
                               struct mdt_properties *prop,
                               const struct mdt_library *mlib,
                               const struct mod_libraries *libs)
{
  if (!prop[is].disulfides) {
    struct mod_sequence *seq = mod_alignment_sequence_get(aln, is);
    struct mod_structure *struc = mod_alignment_structure_get(aln, is);
    prop[is].disulfides = get_disulfides(struc, seq, libs);
  }
  return prop[is].disulfides;
}

/** Require that the tuples have at least min_natom atoms each */
gboolean tuple_require_natom(const struct mdt_library *mlib, int min_natom,
                             GError **err)
{
  int have_natom = mlib->tupclass->natom;
  if (have_natom < min_natom) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "This feature works only for tuples of %d or more "
                "atoms; you have only %d", min_natom, have_natom);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Get a single atom tuple from a structure */
const struct mdt_tuple *property_one_tuple(const struct mod_alignment *aln,
                                           int is, struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int ibnd1, int ia1,
                                           const struct mod_libraries *libs)
{
  const struct mdt_tuple_list *trp;
  trp = property_tuples(aln, is, prop, mlib, libs);

  return &trp[ia1].tuples[ibnd1];
}

const float *property_user(const struct mod_alignment *aln, int is,
                           struct mdt_properties *prop,
                           const struct mdt_library *mlib,
                           const struct mod_libraries *libs, int user_index)
{
  if (!prop[is].user_properties[user_index]) {
    struct mdt_user_property *p;
    p = &g_array_index(mlib->user_properties, struct mdt_user_property,
                       user_index);
    prop[is].user_properties[user_index] = p->get_property(p->data, aln, is,
                                                           mlib, libs);
  }
  return prop[is].user_properties[user_index];
}
