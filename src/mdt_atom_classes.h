/** \file mdt_atom_classes.h  Functions to handle atom classes.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#ifndef __MDT_ATOMCLASS_H
#define __MDT_ATOMCLASS_H

#include <glib.h>
#include "mdt_config.h"
#include "mdt_types.h"

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

/** Make a new atom class list. */
MDTDLLEXPORT
struct mdt_atom_class_list *mdt_atom_class_list_new(int natom);

/** Free an existing atom class list. */
MDTDLLEXPORT
void mdt_atom_class_list_free(struct mdt_atom_class_list *atclass);

/** Read atom class information from a file; return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_atom_classes_read(const gchar *filename,
                               struct mdt_library *mlib, int natom,
                               GError **err);

/** Read hydrogen bond class information from a file; return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_hbond_read(const gchar *filename, struct mdt_library *mlib,
                        GError **err);

/** Read tuple class information from a file; return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_tuple_read(const gchar *filename, struct mdt_library *mlib,
                        GError **err);

/** Set the number of bins and the bin symbols for atom class features */
MDTDLLLOCAL
void update_mdt_feat_atclass(struct mod_mdt_libfeature *feat,
                             const struct mdt_atom_class_list *atclass);

G_END_DECLS

#endif  /* __MDT_ATOMCLASS_H */
