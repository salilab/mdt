/** \file mdt_stereo.c     Functions to determine stereochemistry for
 *                         template structures.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <glib.h>
#include <string.h>
#include "modeller.h"
#include "mdt_atom_classes.h"
#include "mdt_index.h"
#include "mdt_stereo.h"

/** Return True iff the atom index is OK */
gboolean atmdefd(int ia1, const struct mod_coordinates *cd)
{
  return ia1 >= 0 && ia1 < cd->natm;
}

/** Add a bond/angle/dihedral if all the atoms are OK. */
static void add_bond(const struct mod_coordinates *cd,
                     const struct mod_sequence *seq,
                     const struct mdt_atom_type *atype, int ia1, int ir1,
                     int natom, int iclass, GArray *bonds)
{
  struct mdt_bond newbnd;
  int i;

  newbnd.iata[0] = ia1;
  for (i = 1; i < natom; i++) {
    newbnd.iata[i] = mod_residue_find_atom(cd, seq, ir1,
                                           atype->names[i + 1]) - 1;
    if (!atmdefd(newbnd.iata[i], cd)) {
      return;
    }
  }

  /* If all atoms are defined, add the bond to the list */
  newbnd.bndgrp = iclass + 1;
  g_array_append_val(bonds, newbnd);
}

static void get_bondlist(GArray *bonds, const struct mod_structure *struc,
                         const struct mod_sequence *seq,
                         const struct mdt_atom_class_list *atclass,
                         int natom, const struct mod_libraries *libs)
{
  int ia1, *iresatm, *irestyp;
  iresatm = mod_int1_pt(&struc->cd.iresatm);
  irestyp = mod_int1_pt(&seq->irestyp);
  for (ia1 = 0; ia1 < struc->cd.natm; ia1++) {
    int iclass, ir1 = iresatm[ia1] - 1;
    char *resnam = mod_residue_name_from_type(irestyp[ir1], libs);
    for (iclass = 0; iclass < atclass->nclass; iclass++) {
      int itype;
      struct mdt_atom_class *c = &atclass->classes[iclass];
      for (itype = 0; itype < c->ntypes; itype++) {
        struct mdt_atom_type *t = &c->types[itype];

        /* Does the residue type match? */
        if (strcmp(resnam, t->names[0]) == 0 || strcmp(t->names[0], "*") == 0) {
          char *atmnam = mod_coordinates_atmnam_get(&struc->cd, ia1);
          /* Does the lead atom type match? */
          if (strcmp(atmnam, t->names[1]) == 0) {
            add_bond(&struc->cd, seq, t, ia1, ir1, natom, iclass, bonds);
          }
          g_free(atmnam);
        }
      }
    }
    g_free(resnam);
  }
}

/** Get all of one type of bond (bond/angle/dihedral) for a structure. */
struct mdt_bond_list *get_stereo(const struct mod_structure *struc,
                                 const struct mod_sequence *seq,
                                 const struct mdt_atom_class_list *atclass,
                                 int bondtype,
                                 const struct mod_libraries *libs)
{
  struct mdt_bond_list *bondlist = g_malloc(sizeof(struct mdt_bond_list));
  GArray *bonds = g_array_new(FALSE, FALSE, sizeof(struct mdt_bond));

  get_bondlist(bonds, struc, seq, atclass, bondtype + 2, libs);

  bondlist->nbonds = bonds->len;
  bondlist->bonds = (struct mdt_bond *)g_array_free(bonds, FALSE);
  return bondlist;
}
