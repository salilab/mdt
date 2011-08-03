/** \file mdt_residue_bonds.h    Functions to calculate residue bond separation.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#ifndef __MDT_RESIDUE_BONDS_H
#define __MDT_RESIDUE_BONDS_H

#include <glib.h>
#include "mod_types.h"
#include "mdt_config.h"

G_BEGIN_DECLS

/** Number of bonds separating each pair of atoms for a residue type */
struct mdt_residue_bonds {
  /** Mapping from atom names to indices into the distance matrix */
  GHashTable *atom_names;
  /** Number of bonds separating each pair of atoms.
      Note that this is a square matrix; size = g_hash_table_size(atom_names),
      although only the upper triangle is used. */
  int *distance;
};

/** A list of residue bonds */
struct mdt_residue_bond_list {
  /** Number of residue types */
  int nres;
  /** Atom-atom bond distance information for each residue type */
  struct mdt_residue_bonds *bonds;
};

/** Initialize a residue bond list. */
MDTDLLLOCAL
void mdt_residue_bond_list_init(struct mdt_residue_bond_list *bondlist);

/** Free a residue bond list. */
MDTDLLLOCAL
void mdt_residue_bond_list_free(struct mdt_residue_bond_list *bondlist);

struct mdt_library;

/** Fill residue bonds using the bond library. */
MDTDLLLOCAL
void mdt_fill_residue_bonds(struct mdt_residue_bond_list *bondlist,
                            const struct mdt_library *mlib,
                            struct mod_libraries *libs);

G_END_DECLS

#endif  /* __MDT_RESIDUE_BONDS_H */
