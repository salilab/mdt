/** \file mdt_types.h      Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_TYPES_H
#define __MDT_TYPES_H

#include "mdt_config.h"
#include "mod_mdt_type.h"

#include <glib.h>

G_BEGIN_DECLS

/** Atom type */
struct mdt_atom_type {
  /** Atom names */
  char **names;
};

/** Atom class, which can contain multiple atom types */
struct mdt_atom_class {
  /** Number of atom types */
  int ntypes;
  /** Atom types */
  struct mdt_atom_type *types;
  /** Class name */
  char *name;
  /** Hydrogen bond properties (donor/acceptor valency, and charge) */
  float hb_property[3];
};

/** List of atom classes */
struct mdt_atom_class_list {
  /** Number of atoms in each atom type (1=atom, 2=bond, 3=angle, 4=dihedral) */
  int natom;
  /** Number of classes */
  int nclass;
  /** Classes */
  struct mdt_atom_class *classes;
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
};

/** User-defined feature types */
typedef enum {
  MDT_FEATURE_NONE = 0,
  MDT_FEATURE_PROTEIN,
  MDT_FEATURE_RESIDUE,
  MDT_FEATURE_ATOM
} mdt_feature_type;

struct mdt_properties;

typedef int (*mdt_cb_feature_protein)(const struct mod_alignment *aln,
                                      int protein,
                                      struct mdt_properties *prop, void *data,
                                      const struct mod_mdt_libfeature *feat,
                                      const struct mod_libraries *libs,
                                      GError **err);


typedef int (*mdt_cb_feature_residue)(const struct mod_alignment *aln,
                                      int protein, int residue,
                                      struct mdt_properties *prop, void *data,
                                      const struct mod_mdt_libfeature *feat,
                                      const struct mod_libraries *libs,
                                      GError **err);

typedef int (*mdt_cb_feature_atom)(const struct mod_alignment *aln,
                                   int protein, int atom,
                                   struct mdt_properties *prop, void *data,
                                   const struct mod_mdt_libfeature *feat,
                                   const struct mod_libraries *libs,
                                   GError **err);

/** User-defined protein feature */
struct mdt_feature_protein {
  int protein;
  mdt_cb_feature_protein getbin;
};

/** User-defined residue feature */
struct mdt_feature_residue {
  int protein;
  int delta;
  gboolean pos2;
  mdt_cb_feature_residue getbin;
};

/** User-defined atom feature */
struct mdt_feature_atom {
  gboolean pos2;
  mdt_cb_feature_atom getbin;
};

/** User-defined feature */
struct mdt_feature {
  mdt_feature_type type;
  union {
    struct mdt_feature_protein protein;
    struct mdt_feature_residue residue;
    struct mdt_feature_atom atom;
  } u;
  void *data;
};

/** Library of feature data used by MDTs */
struct mdt_library {
  /** Base Modeller type */
  struct mod_mdt_library base;
  /** Deltas for some feature types ('residue + delta i') */
  int deltai, deltaj;
  /** TRUE if deltas refer to align. positions, or FALSE if residue positions */
  gboolean deltai_ali, deltaj_ali;
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
};

/** Make a new mdt structure */
MDTDLLEXPORT
struct mdt *mdt_new(mod_mdt_bin_type bin_type);

/** Free an mdt structure */
MDTDLLEXPORT
void mdt_free(struct mdt *mdt);

/** Make a new mdt_library structure */
MDTDLLEXPORT
struct mdt_library *mdt_library_new(void);

/** Free an mdt_library structure */
MDTDLLEXPORT
void mdt_library_free(struct mdt_library *mlib);

G_END_DECLS

#endif  /* __MDT_TYPES_H */
