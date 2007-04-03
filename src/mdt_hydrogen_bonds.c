/** \file mdt_hydrogen_bonds.c    Functions for handling hydrogen bonds
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <math.h>
#include <glib.h>
#include "mdt_types.h"
#include "mdt_hydrogen_bonds.h"
#include "modeller.h"

static float dist1(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float xdiff, ydiff, zdiff;
  xdiff = x1 - x2;
  ydiff = y1 - y2;
  zdiff = z1 - z2;
  return sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
}

static float get_ndon_acc(const struct mdt_atom_class_list *atclass, int iat,
                          gboolean acceptor)
{
  int i, j;
  for (i = 0; i < atclass->nclass; i++) {
    struct mdt_atom_class *c = &atclass->classes[i];
    for (j = 0; j < c->ntypes; j++) {
      struct mdt_atom_type *t = &c->types[j];
      iat--;
      if (iat < 0) {
        return (acceptor ? ABS(t->hb_acceptor) : ABS(t->hb_donor));
      }
    }
  }
  return 0;
}

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
int numb_hda(int ia, const int hb_iattyp[], const struct coordinates *cd,
             const struct mdt_atom_class_list *atclass, float hbond_cutoff,
             gboolean acceptor, int nbins)
{
  float *x, *y, *z;
  int i, num = 0;
  x = f_float1_pt(&cd->x);
  y = f_float1_pt(&cd->y);
  z = f_float1_pt(&cd->z);
  ia--;
  for (i = 0; i < cd->natm; i++) {
    if (i != ia) {
      int iat = hb_iattyp[i] - 1;
      if (iat >= 0 && get_ndon_acc(atclass, iat, acceptor) > 0.) {
        float d = dist1(x[i], y[i], z[i], x[ia], y[ia], z[ia]);
        if (d > 2.5 && d < hbond_cutoff) {
          num++;
        }
      }
    }
  }
  
  return MIN(num, nbins);
}
