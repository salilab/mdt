/** \file mdt_types.h      Functions to handle MDT types.
 *
 *             Part of MODELLER, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_TYPES_H
#define __MDT_TYPES_H

#include <glib.h>
#include "mod_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Atom type */
struct mdt_atom_type {
  /** Atom names */
  char **names;
  /* Hydrogen bond donor valency */
  float hb_donor;
  /* Hydrogen bond acceptor valency */
  float hb_acceptor;
  /* Hydrogen bond charge */
  float hb_charge;
};

/** Atom class, which can contain multiple atom types */
struct mdt_atom_class {
  /** Number of atom types */
  int ntypes;
  /** Atom types */
  struct mdt_atom_type *types;
  /** Class name */
  char *name;
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

/** Library of feature data used by MDTs */
struct mdt_library {
  /** Base Modeller type */
  struct mod_mdt_library base;
  /** Deltas for some feature types ('residue + delta i') */
  int deltai, deltaj;
  /** TRUE if deltas refer to align. positions, or FALSE if residue positions */
  gboolean deltai_ali, deltaj_ali;
  /** Atom, bond, angle, dihedral classes */
  struct mdt_atom_class_list *atclass[4];
};

/** Make a new mdt_library structure */
struct mdt_library *mdt_library_new(void);

/** Free an mdt_library structure */
void mdt_library_free(struct mdt_library *mlib);

G_END_DECLS

#endif  /* __MDT_TYPES_H */
