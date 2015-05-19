/** \file mdt_types.h      Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2015 Andrej Sali
 */

#ifndef __MDT_TYPES_H
#define __MDT_TYPES_H

#include "mdt_config.h"
#include "mod_mdt_type.h"
#include "mod_types.h"
#include "mdt_residue_bonds.h"
#include "mdt_hdf5.h"

#include <glib.h>

G_BEGIN_DECLS

/** Modeller dihedral angle types.
    These correspond to types hard-coded into Modeller and defined in
    Modeller's modlib/resdih.lib file, so should not be changed. */
typedef enum {
  MDT_DIHEDRAL_ALPHA = 0,
  MDT_DIHEDRAL_PHI,
  MDT_DIHEDRAL_PSI,
  MDT_DIHEDRAL_OMEGA,
  MDT_DIHEDRAL_CHI1,
  MDT_DIHEDRAL_CHI2,
  MDT_DIHEDRAL_CHI3,
  MDT_DIHEDRAL_CHI4,
  MDT_DIHEDRAL_CHI5
} mdt_dihedral_type;

struct mdt_scan_parameters {
  /** True only after at least one scan has been performed */
  gboolean scan_called;
  int residue_span_range[4];
  int chain_span_range[4];
  int bond_span_range[2];
  gboolean disulfide, exclude_bonds, exclude_angles, exclude_dihedrals,
           sympairs, symtriples;
  float distngh;
  int surftyp, accessibility_type;
};

/** A single MDT */
struct mdt {
  /** Base Modeller type */
  struct mod_mdt base;
  /** Number of alignments used to get this MDT */
  int nalns;
  /** Number of prots/prot pairs in analyzed alignments */
  int n_protein_pairs;
  /** Number of proteins in analyzed alignments */
  int n_proteins;
  /** Sample points in the current MDT */
  double sample_size;
  /** TRUE if the MDT is a pdf, FALSE if frequencies */
  gboolean pdf;
  /** TRUE if all features are symmetric */
  gboolean symmetric;
  /** Scan type */
  int scantype;
  /** Callbacks for writing information used by this table to HDF5 files */
  GHashTable *write_lib_funcs;
  /** Parameters used in last scan */
  struct mdt_scan_parameters scan_params;
};

struct mdt_atom_class_list;

/** Function to get a property for a structure; it should return a malloc'd
    float array, which MDT will free when it is no longer needed. On error,
    it should return NULL. */
typedef float * (*mdt_cb_get_property)(gpointer data,
                                       const struct mod_alignment *aln, int is,
                                       const struct mdt_library *mlib,
                                       const struct mod_libraries *libs,
                                       GError **err);

/** Callbacks to populate a user-defined property */
struct mdt_user_property {
  mdt_cb_get_property get_property;
  gpointer data;
  GDestroyNotify freefunc;
};

/** Library of feature data used by MDTs */
struct mdt_library {
  /** Base Modeller type */
  struct mod_mdt_library base;
  /** Modeller libraries */
  struct mod_libraries *libs;
  /** Scripting language object */
  gpointer scriptobj;
  /** Whether to treat disulfides and termini specially for atom types */
  gboolean special_atoms;
  /** Cutoff distance for hydrogen bonds */
  float hbond_cutoff;
  /** User-defined features */
  GArray *features;
  /** Atom, bond, angle, dihedral classes */
  struct mdt_atom_class_list *atclass[4];
  /** Hydrogen bond classes */
  struct mdt_atom_class_list *hbond;
  /** Tuple classes */
  struct mdt_atom_class_list *tupclass;
  /** Set TRUE if at least one feature has been added. */
  gboolean feature_added;
  /** Atom names for distance calculation */
  char *distance_atoms[2];
  /** For each residue, number of bonds separating each pair of atoms */
  struct mdt_residue_bond_list residue_bond_list;
  /** User-defined properties (struct mdt_user_property) */
  GArray *user_properties;
};

typedef gboolean (*mdt_cb_write_lib)(hid_t loc_id,
                                     const struct mdt_library *mlib);

/** Make a new mdt structure */
MDTDLLEXPORT
struct mdt *mdt_new(mod_mdt_bin_type bin_type);

/** Register a callback function to write data used by this table to HDF5 */
MDTDLLEXPORT
void mdt_set_write_lib_callback(struct mdt *mdt, mdt_cb_write_lib writelibfunc);

/** Free an mdt structure */
MDTDLLEXPORT
void mdt_free(struct mdt *mdt);

/** Make a new mdt_library structure */
MDTDLLEXPORT
struct mdt_library *mdt_library_new(struct mod_libraries *libs,
                                    gpointer scriptobj);

/** Free an mdt_library structure */
MDTDLLEXPORT
void mdt_library_free(struct mdt_library *mlib);

/** Register a user-defined property. Return the ID of that property. */
MDTDLLEXPORT
int mdt_library_add_user_property(struct mdt_library *mlib,
                                  mdt_cb_get_property get_property,
                                  gpointer data, GDestroyNotify freefunc);

G_END_DECLS

#endif  /* __MDT_TYPES_H */
