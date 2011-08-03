/** \file mdt_residue_bonds.c    Functions to calculate residue bond separation.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <glib.h>
#include "mdt_types.h"
#include "mdt_residue_bonds.h"
#include "mdt_atom_classes.h"
#include "modeller.h"

/** Initialize a residue bond list. */
void mdt_residue_bond_list_init(struct mdt_residue_bond_list *bondlist)
{
  bondlist->nres = 0;
  bondlist->bonds = NULL;
}

/** Free a residue bond list. */
void mdt_residue_bond_list_free(struct mdt_residue_bond_list *bondlist)
{
  int i;
  for (i = 0; i < bondlist->nres; ++i) {
    struct mdt_residue_bonds *resbonds = &bondlist->bonds[i];
    if (resbonds->atom_names) {
      g_hash_table_destroy(resbonds->atom_names);
      g_free(resbonds->distance);
    }
  }
}

/* Do initial setup of the residue bond list. */
static GHashTable *make_residue_types(struct mdt_residue_bond_list *bondlist,
                                      const int nres,
                                      struct mod_libraries *libs)
{
  GHashTable *res_hash;
  int i;

  res_hash = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  for (i = 0; i < nres; ++i) {
    g_hash_table_insert(res_hash, mod_residue_name_from_type(i + 1, libs),
                        GINT_TO_POINTER(i + 1));
  }

  bondlist->nres = nres;
  bondlist->bonds = g_malloc(sizeof(struct mdt_residue_bonds) * nres);
  for (i = 0; i < nres; ++i) {
    bondlist->bonds[i].distance = NULL;
    bondlist->bonds[i].atom_names = NULL;
  }
  return res_hash;
}

/* Is the bond between the two named atoms confined to a single residue,
   as opposed to spanning residues? */
static gboolean internal_bond(const char *at1, const char *at2)
{
  return (at1[0] != '+' && at1[0] != '-' && at2[0] != '+' && at2[0] != '-');
}

/* Add the named atom to the list of atoms in this residue */
static void add_residue_atom(struct mdt_residue_bonds *resbond, char *atom)
{
  if (!resbond->atom_names) {
    /* Note: we don't own the atom string pointers (the atom class
       struct does) so don't free them when we're done with the hash. */
    resbond->atom_names = g_hash_table_new(g_str_hash, g_str_equal);
    /* Add backbone N and C atoms at start of table for speed */
    g_hash_table_insert(resbond->atom_names, "N", GINT_TO_POINTER(0));
    g_hash_table_insert(resbond->atom_names, "C", GINT_TO_POINTER(1));
  }

  if (!g_hash_table_lookup_extended(resbond->atom_names, atom, NULL, NULL)) {
    int index = g_hash_table_size(resbond->atom_names);
    g_hash_table_insert(resbond->atom_names, atom, GINT_TO_POINTER(index));
  }
}

/* Populate the list of atoms for every residue type */
static void make_atom_names(struct mdt_residue_bonds *resbonds,
                            struct mdt_atom_class_list *bondcl,
                            GHashTable *res_hash)
{
  int icls, ityp;
  for (icls = 0; icls < bondcl->nclass; ++icls) {
    struct mdt_atom_class *cls = &bondcl->classes[icls];
    for (ityp = 0; ityp < cls->ntypes; ++ityp) {
      struct mdt_atom_type *typ = &cls->types[ityp];
      int res = GPOINTER_TO_INT(g_hash_table_lookup(res_hash, typ->names[0]));
      if (res && internal_bond(typ->names[1], typ->names[2])) {
        add_residue_atom(&resbonds[res-1], typ->names[1]);
        add_residue_atom(&resbonds[res-1], typ->names[2]);
      }
    }
  }
}

/* Set the bond distance between the two atom types in the residue */
static void set_distance(struct mdt_residue_bonds *resbond, int at1, int at2,
                         int distance)
{
  int natom = g_hash_table_size(resbond->atom_names);
  int index = MIN(at1, at2) * natom + MAX(at1, at2);
  resbond->distance[index] = distance;
}

/* Get the bond distance between the two atom types in the residue */
static int get_distance(struct mdt_residue_bonds *resbond, int at1, int at2)
{
  int natom = g_hash_table_size(resbond->atom_names);
  int index = MIN(at1, at2) * natom + MAX(at1, at2);
  return resbond->distance[index];
}

/* Follow bonds recursively to find distances to all atoms from the given
   start atom */
static void follow_bond(struct mdt_residue_bonds *resbond, int natom,
                        int start_atom, int endpoint1, int endpoint2,
                        int distance)
{
  int other_atom;
  /* Protect against infinite loops if cycles are present in the structure */
  if (distance > natom) {
    return;
  }
  for (other_atom = 0; other_atom < natom; ++other_atom) {
    /* Don't double back on ourselves */
    if (other_atom != start_atom && other_atom != endpoint1
        && other_atom != endpoint2) {
      int dist = get_distance(resbond, endpoint2, other_atom);
      if (dist == 1) {
        int old_distance = get_distance(resbond, start_atom, other_atom);
        /* In case of cycles, take the shortest route */
        if (distance + 1 <= old_distance) {
          set_distance(resbond, start_atom, other_atom, distance + 1);
          follow_bond(resbond, natom, start_atom, endpoint2, other_atom,
                      distance + 1);
        }
      }
    }
  }
}

/* Get bond distances between all pairs of atoms in every residue type */
static void fill_distances(struct mdt_residue_bond_list *bondlist)
{
  int ires;

  for (ires = 0; ires < bondlist->nres; ++ires) {
    struct mdt_residue_bonds *resbond = &bondlist->bonds[ires];
    int natom = g_hash_table_size(resbond->atom_names);
    int start_atom;
    for (start_atom = 0; start_atom < natom; ++start_atom) {
      int other_atom;
      for (other_atom = 0; other_atom < natom; ++other_atom) {
        if (other_atom != start_atom) {
          int dist = get_distance(resbond, start_atom, other_atom);
          if (dist == 1) {
            follow_bond(resbond, natom, start_atom, -1, other_atom, 1);
          }
        }
      }
    }
  }
}

/* Initialize the atom-atom bond distance matrix */
static void init_distances(struct mdt_residue_bond_list *bondlist)
{
  int i, j, k;
  for (i = 0; i < bondlist->nres; ++i) {
    struct mdt_residue_bonds *resbonds = &bondlist->bonds[i];
    if (resbonds->atom_names) {
      int ind = 0;
      int natom = g_hash_table_size(resbonds->atom_names);
      resbonds->distance = g_malloc(sizeof(int) * natom * natom);
      for (j = 0; j < natom; ++j) {
        for (k = 0; k < natom; ++k) {
          resbonds->distance[ind++] = (k == j ? 0 : 10000);
        }
      }
    }
  }
}

/* Transfer bond information to the atom-atom distance matrix */
static void add_bonds(struct mdt_residue_bonds *resbonds,
                      struct mdt_atom_class_list *bondcl,
                      GHashTable *res_hash)
{
  int icls, ityp;

  for (icls = 0; icls < bondcl->nclass; ++icls) {
    struct mdt_atom_class *cls = &bondcl->classes[icls];
    for (ityp = 0; ityp < cls->ntypes; ++ityp) {
      struct mdt_atom_type *typ = &cls->types[ityp];
      int res = GPOINTER_TO_INT(g_hash_table_lookup(res_hash, typ->names[0]));
      if (res && internal_bond(typ->names[1], typ->names[2])) {
        int at1, at2;
        at1 = GPOINTER_TO_INT(g_hash_table_lookup(resbonds[res-1].atom_names,
                                                  typ->names[1]));
        at2 = GPOINTER_TO_INT(g_hash_table_lookup(resbonds[res-1].atom_names,
                                                  typ->names[2]));
        set_distance(&resbonds[res-1], at1, at2, 1);
      }
    }
  }
}

/** Fill residue bonds using the bond library. */
void mdt_fill_residue_bonds(struct mdt_residue_bond_list *bondlist,
                            const struct mdt_library *mlib,
                            struct mod_libraries *libs)
{
  /* Consider only standard amino acids */
  const int n_std_restyp = 20;

  GHashTable *res_hash;

  /* Do nothing if the information is already populated */
  if (bondlist->nres > 0) {
    return;
  }

  res_hash = make_residue_types(bondlist, n_std_restyp, libs);
  make_atom_names(bondlist->bonds, mlib->atclass[1], res_hash);
  init_distances(bondlist);
  add_bonds(bondlist->bonds, mlib->atclass[1], res_hash);
  g_hash_table_destroy(res_hash);
  fill_distances(bondlist);
}
