/** \file mdt_tuples.c     Functions to build lists of atom tuples.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
 */

#include <string.h>
#include <glib.h>
#include "modeller.h"
#include "mdt_types.h"
#include "mdt_index.h"
#include "mdt_stereo.h"
#include "mdt_tuples.h"
#include "mdt_atom_classes.h"

static int find_atom(const struct mod_coordinates *cd,
                     const struct mod_sequence *seq, int ir1,
                     const char *name, const int *irestyp,
                     const struct mod_libraries *libs)
{
  if (name[0] == '(') {
    char *right = strchr(name, ')');
    if (right && right[1] != '\0') {
      int atind = mod_residue_find_atom(cd, seq, ir1, &right[1]) - 1;
      if (atind >= 0) {
        int cmp;
        char *resnam;
        if (right[1] == '+') {
          ir1++;
        } else if (right[1] == '-') {
          ir1--;
        }
        resnam = mod_residue_name_from_type(irestyp[ir1], libs);
        cmp = strncmp(&name[1], resnam, right - name - 1);
        g_free(resnam);
        if (cmp != 0) {
          return -1;
        }
      }
      return atind;
    }
  }
  return mod_residue_find_atom(cd, seq, ir1, name) - 1;
}

static void add_tuple(const struct mod_coordinates *cd,
                      const struct mod_sequence *seq,
                      const struct mdt_atom_type *atype, int ia1, int ir1,
                      int iclass, int natom, const int *irestyp,
                      const struct mod_libraries *libs, GArray *tlist)
{
  int ia2, ia3 = 0;
  struct mdt_tuple newtup;

  /* 2nd atom of the tuple: */
  ia2 = find_atom(cd, seq, ir1, atype->names[2], irestyp, libs);
  if (ia1 == ia2 || !atmdefd(ia2, cd)) {
    return;
  }
  if (natom > 2) {
    /* 3rd atom of the tuple (if applicable): */
    ia3 = find_atom(cd, seq, ir1, atype->names[3], irestyp, libs);
    if (ia1 == ia3 || ia2 == ia3 || !atmdefd(ia3, cd)) {
      return;
    }
  }

  newtup.iata[0] = ia2;
  newtup.iata[1] = ia3;
  newtup.tupclass = iclass + 1;
  g_array_append_val(tlist, newtup);
}

static void get_tuples(int iatm, const struct mod_structure *struc,
                       const struct mod_sequence *seq,
                       const struct mdt_atom_class_list *atclass,
                       GArray *tlist, const struct mod_libraries *libs)
{
  int iclass, *iresatm, *irestyp, ir1;
  char *resnam;

  iresatm = mod_int1_pt(&struc->cd.iresatm);
  irestyp = mod_int1_pt(&seq->irestyp);
  ir1 = iresatm[iatm] - 1;
  resnam = mod_residue_name_from_type(irestyp[ir1], libs);

  for (iclass = 0; iclass < atclass->nclass; iclass++) {
    int itype;
    struct mdt_atom_class *c = &atclass->classes[iclass];
    for (itype = 0; itype < c->ntypes; itype++) {
      struct mdt_atom_type *t = &c->types[itype];

      /* Does the residue type match? */
      if (strcmp(resnam, t->names[0]) == 0 || strcmp(t->names[0], "*") == 0) {
        char *atmnam = mod_coordinates_atmnam_get(&struc->cd, iatm);
        /* Does the lead atom type match? */
        if (strcmp(atmnam, t->names[1]) == 0) {
          add_tuple(&struc->cd, seq, t, iatm, ir1, iclass, atclass->natom,
                    irestyp, libs, tlist);
        }
        g_free(atmnam);
      }
    }
  }
  g_free(resnam);
}


/** Get all tuples for a structure. */
struct mdt_tuple_list *tupclass(const struct mod_structure *struc,
                                const struct mod_sequence *seq,
                                const struct mdt_atom_class_list *atclass,
                                const struct mod_libraries *libs)
{
  int i;
  struct mdt_tuple_list *tuples;
  GArray *tlist;
  tuples = g_malloc(sizeof(struct mdt_tuple_list) * struc->cd.natm);
  tlist = g_array_new(FALSE, FALSE, sizeof(struct mdt_tuple));
  for (i = 0; i < struc->cd.natm; i++) {
    if (atmdefd(i, &struc->cd)) {
      get_tuples(i, struc, seq, atclass, tlist, libs);
    }
    if (tlist->len > 0) {
      tuples[i].ntuples = tlist->len;
      tuples[i].tuples = (struct mdt_tuple *)g_array_free(tlist, FALSE);
      tlist = g_array_new(FALSE, FALSE, sizeof(struct mdt_tuple));
    } else {
      tuples[i].ntuples = 0;
      tuples[i].tuples = NULL;
    }
  }
  g_array_free(tlist, TRUE);
  return tuples;
}
