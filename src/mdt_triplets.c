/** \file mdt_triplets.c   Functions to build lists of atom triplets.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <string.h>
#include <glib.h>
#include "modeller.h"
#include "mdt_types.h"
#include "mdt_index.h"
#include "mdt_stereo.h"
#include "mdt_triplets.h"

static void add_triplet(const struct coordinates *cd,
                        const struct sequence *seq,
                        const struct mdt_atom_type *atype, int ia1, int ir1,
                        int iclass, int natom, GArray *tlist)
{
  int ia2, ia3 = 0;
  struct mdt_triplet newtrp;

  /* 2nd atom of the triplet: */
  ia2 = mod_residue_find_atom(cd, seq, ir1, atype->names[2]) - 1;
  if (ia1 == ia2 || !atmdefd(ia2, cd)) {
    return;
  }
  if (natom > 2) {
    /* 3rd atom of the triplet: */
    ia3 = mod_residue_find_atom(cd, seq, ir1, atype->names[3]) - 1;
    if (ia1 == ia3 || ia2 == ia3 || !atmdefd(ia3, cd)) {
      return;
    }
  }

  newtrp.iata[0] = ia2;
  newtrp.iata[1] = ia3;
  newtrp.trpclass = iclass + 1;
  g_array_append_val(tlist, newtrp);
}

static void get_triplets(int iatm, const struct structure *struc,
                         const struct sequence *seq,
                         const struct mdt_atom_class_list *atclass,
                         GArray *tlist, const struct libraries *libs)
{
  int iclass, *iresatm, *irestyp, ir1;
  char *resnam;

  iresatm = f_int1_pt(&struc->cd.iresatm);
  irestyp = f_int1_pt(&seq->irestyp);
  ir1 = iresatm[iatm] - 1;
  resnam = residue_name_from_type(irestyp[ir1], libs);

  for (iclass = 0; iclass < atclass->nclass; iclass++) {
    int itype;
    struct mdt_atom_class *c = &atclass->classes[iclass];
    for (itype = 0; itype < c->ntypes; itype++) {
      struct mdt_atom_type *t = &c->types[itype];

      /* Does the residue type match? */
      if (strcmp(resnam, t->names[0]) == 0 || strcmp(t->names[0], "*") == 0) {
        char *atmnam = get_coord_atmnam(&struc->cd, iatm);
        /* Does the lead atom type match? */
        if (strcmp(atmnam, t->names[1]) == 0) {
          add_triplet(&struc->cd, seq, t, iatm, ir1, iclass, atclass->natom,
                      tlist);
        }
        g_free(atmnam);
      }
    }
  }
  g_free(resnam);
}


/** Get all triplets for a structure. */
struct mdt_triplet_list *trpclass(const struct structure *struc,
                                  const struct sequence *seq,
                                  const struct mdt_atom_class_list *atclass,
                                  const struct libraries *libs)
{
  int i;
  struct mdt_triplet_list *triplets;
  GArray *tlist;
  triplets = g_malloc(sizeof(struct mdt_triplet_list) * struc->cd.natm);
  tlist = g_array_new(FALSE, FALSE, sizeof(struct mdt_triplet));
  for (i = 0; i < struc->cd.natm; i++) {
    if (atmdefd(i, &struc->cd)) {
      get_triplets(i, struc, seq, atclass, tlist, libs);
    }
    if (tlist->len > 0) {
      triplets[i].ntriplets = tlist->len;
      triplets[i].triplets = (struct mdt_triplet *)g_array_free(tlist, FALSE);
      tlist = g_array_new(FALSE, FALSE, sizeof(struct mdt_triplet));
    } else {
      triplets[i].ntriplets = 0;
      triplets[i].triplets = NULL;
    }
  }
  g_array_free(tlist, TRUE);
  return triplets;
}
