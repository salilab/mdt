/** \file mdt_types.h      Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#ifndef __MDT_TYPES_H
#define __MDT_TYPES_H

#include "mdt_config.h"
#include "mod_mdt_type.h"
#include "mod_types.h"
#include "mdt_residue_bonds.h"

#include <glib.h>

G_BEGIN_DECLS

/** Modeller dihedral angle types */
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

struct mdt_atom_class_list;

/** Library of feature data used by MDTs */
struct mdt_library {
  /** Base Modeller type */
  struct mod_mdt_library base;
  /** Modeller libraries */
  const struct mod_libraries *libs;
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
};

/** Make a new mdt structure */
MDTDLLEXPORT
struct mdt *mdt_new(mod_mdt_bin_type bin_type);

/** Free an mdt structure */
MDTDLLEXPORT
void mdt_free(struct mdt *mdt);

/** Make a new mdt_library structure */
MDTDLLEXPORT
struct mdt_library *mdt_library_new(const struct mod_libraries *libs);

/** Free an mdt_library structure */
MDTDLLEXPORT
void mdt_library_free(struct mdt_library *mlib);

G_END_DECLS

#endif  /* __MDT_TYPES_H */
