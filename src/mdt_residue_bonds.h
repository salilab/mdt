/** \file mdt_residue_bonds.h    Functions to calculate residue bond separation.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#ifndef __MDT_RESIDUE_BONDS_H
#define __MDT_RESIDUE_BONDS_H

#include <glib.h>
#include "mod_types.h"
#include "mdt_config.h"


G_BEGIN_DECLS

struct mdt_disulfide_list;

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
                            const struct mod_libraries *libs);

/** Assign atom types to a structure.
    The array of atom types is returned. It is the caller's responsibility
    to free it when it is no longer needed. */
MDTDLLLOCAL
int *mdt_residue_bonds_assign_atom_types(const struct mod_structure *struc,
                        const struct mod_sequence *seq,
                        const struct mdt_residue_bond_list *bondlist,
                        const struct mod_libraries *libs);

/** Get the number of bonds separating two atoms in a structure.
    -1 is returned if the atoms are not connected. */
MDTDLLLOCAL
int mdt_get_bond_separation(const struct mod_structure *struc,
                            const struct mod_sequence *seq,
                            int atom1, int atom2, const int *attyp,
                            const struct mdt_residue_bond_list *bondlist,
                            const struct mdt_disulfide_list *disulfides);

/** Get the number of bonds separating two atoms in the same chain.
    -1 is returned if the atoms are not connected.
 */
int mdt_get_bond_separation_same_chain(int atom1, int atom2, int res1,
                            int res2, const struct mod_sequence *seq,
                            const int *attyp,
                            const struct mdt_residue_bond_list *bondlist,
                            const struct mdt_disulfide_list *disulfides);


G_END_DECLS

#endif  /* __MDT_RESIDUE_BONDS_H */
