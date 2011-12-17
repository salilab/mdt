/** \file mdt_disulfides.c    Functions to handle disulfide bridges.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <glib.h>
#include "mdt_types.h"
#include "mdt_disulfides.h"
#include "mdt_property.h"
#include "modeller.h"

/** Get list of disulfides. */
void get_disulfides(const struct mod_structure *struc,
                    const struct mod_sequence *seq,
                    const struct mod_libraries *libs,
                    struct mdt_properties *prop, int is)
{
  const char *routine = "get_disulfides";
  int ires, i, j;
  GArray *sg_atoms;
  int *issa = NULL, *issap = NULL, *iresatm;
  int numofss = 0;
  float *x, *y, *z;
  int cystyp = mod_residue_type_from_name("CYS", libs);

  if (cystyp == 0) {
    mod_logwarning(routine, "No CYS residue type; disulfides not detected");
    return;
  }

  sg_atoms = g_array_new(FALSE, FALSE, sizeof(int));
  for (ires = 0; ires < seq->nres; ++ires) {
    int iat = mod_residue_find_atom(&struc->cd, seq, ires, "SG") - 1;
    if (iat >= 0) {
      g_array_append_val(sg_atoms, iat);
    }
  }

  /* Find the ss bond atoms */
  x = mod_float1_pt(&struc->cd.x);
  y = mod_float1_pt(&struc->cd.y);
  z = mod_float1_pt(&struc->cd.z);
  for (i = 0; i < sg_atoms->len; ++i) {
    int iat = g_array_index(sg_atoms, int, i);
    for (j = i + 1; j < sg_atoms->len; ++j) {
      int jat = g_array_index(sg_atoms, int, j);
      float xd, yd, zd, dist2;
      xd = x[iat] - x[jat];
      yd = y[iat] - y[jat];
      zd = z[iat] - z[jat];
      dist2 = xd * xd + yd * yd + zd * zd;
      if (dist2 < 2.5 * 2.5) {
        numofss++;
        issa = g_realloc(issa, sizeof(int) * numofss);
        issa[numofss - 1] = iat;
        issap=g_realloc(issap, sizeof(int) * numofss);
        issap[numofss - 1] = jat;
      }
    }
  }
  prop[is].issa=issa;
  prop[is].issap=issap;
  prop[is].numofss=numofss;
  iresatm = mod_int1_pt(&struc->cd.iresatm);

  if (numofss>0){
    mod_lognote("%d s-s bond",numofss);
    prop[is].issr=g_malloc(sizeof(int)*numofss);
    prop[is].issrp=g_malloc(sizeof(int)*numofss);
    for(i = 0; i < numofss; i++){
      prop[is].issr[i]=iresatm[prop[is].issa[i]]-1;
      prop[is].issrp[i]=iresatm[prop[is].issap[i]]-1;
      mod_lognote("Residue %s ATOM1 %s %d %d, ATOM2, %s %d %d",
                  mod_residue_name_from_type(2,libs),
                  mod_coordinates_atmnam_get(&struc->cd,prop[is].issa[i]),
                  prop[is].issr[i], prop[is].issa[i],
                  mod_coordinates_atmnam_get(&struc->cd,prop[is].issa[i]),
                  prop[is].issrp[i], prop[is].issap[i]);
    }
  }

  g_array_free(sg_atoms, TRUE);
}
