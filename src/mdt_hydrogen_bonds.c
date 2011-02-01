/** \file mdt_hydrogen_bonds.c    Functions for handling hydrogen bonds
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <math.h>
#include <glib.h>
#include "mdt_types.h"
#include "mdt_atom_classes.h"
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

static float get_hbprop(const struct mdt_atom_class_list *atclass, int iat,
                        int hbprop_type)
{
  struct mdt_atom_class *c = &atclass->classes[iat];
  return ABS(c->hb_property[hbprop_type]);
}

/** Return the number of H-bonds with a given atom, ia. */
int numb_hda(int ia, const int hb_iattyp[], const struct mod_coordinates *cd,
             const struct mdt_atom_class_list *atclass, float hbond_cutoff,
             int hbprop_type)
{
  float *x, *y, *z;
  int i, num = 0;
  x = mod_float1_pt(&cd->x);
  y = mod_float1_pt(&cd->y);
  z = mod_float1_pt(&cd->z);
  for (i = 0; i < cd->natm; i++) {
    if (i != ia) {
      int iat = hb_iattyp[i] - 1;
      if (iat >= 0 && get_hbprop(atclass, iat, hbprop_type) > 0.) {
        float d = dist1(x[i], y[i], z[i], x[ia], y[ia], z[ia]);
        if (d > 2.5 && d < hbond_cutoff) {
          num++;
        }
      }
    }
  }

  return num;
}

/** Calculate H-bond protein satisfaction for the whole protein. */
float hb_satisfaction(const struct mod_coordinates *cd, const int hb_iattyp[],
                      const struct mdt_atom_class_list *atclass,
                      float hbond_cutoff)
{
  float satis = 0.;
  int ia;

  for (ia = 0; ia < cd->natm; ia++) {
    int iprop;
    int iat = hb_iattyp[ia] - 1;
    /* Loop over donors and acceptors */
    for (iprop = 0; iprop <= 1; iprop++) {
      float hbprop;
      int nhda;
      hbprop = atclass->classes[iat].hb_property[iprop];
      nhda = numb_hda(ia, hb_iattyp, cd, atclass, hbond_cutoff, iprop);
      satis += MAX(0., hbprop - nhda);
    }
  }

  return (cd->natm > 0 ? satis / cd->natm : 0.);
}
