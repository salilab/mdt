/** \file mdt_disulfides.c    Functions to handle disulfide bridges.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <glib.h>
#include "mdt_types.h"
#include "mdt_disulfides.h"
#include "mdt_property.h"
#include "modeller.h"

/** Find all SG atoms in CYS residues */
static GArray *find_sg_atoms(const struct mod_structure *struc,
                             const struct mod_sequence *seq,
                             const struct mod_libraries *libs)
{
  int ires;
  int cystyp = mod_residue_type_from_name("CYS", libs);
  GArray *sg_atoms = g_array_new(FALSE, FALSE, sizeof(int));
  int *irestyp = mod_int1_pt(&seq->irestyp);
  for (ires = 0; ires < seq->nres; ++ires) {
    if (irestyp[ires] == cystyp) {
      int iat = mod_residue_find_atom(&struc->cd, seq, ires, "SG") - 1;
      if (iat >= 0) {
        g_array_append_val(sg_atoms, iat);
      }
    }
  }
  return sg_atoms;
}

/** Get list of disulfides. */
struct mdt_disulfide_list *get_disulfides(const struct mod_structure *struc,
                                          const struct mod_sequence *seq,
                                          const struct mod_libraries *libs)
{
  struct mdt_disulfide_list *disulfides;
  int i, j;
  GArray *sg_atoms, *ss;
  int *iresatm;
  float *x, *y, *z;

  sg_atoms = find_sg_atoms(struc, seq, libs);

  /* Find the ss bond atoms */
  ss = g_array_new(FALSE, FALSE, sizeof(struct mdt_disulfide));
  x = mod_float1_pt(&struc->cd.x);
  y = mod_float1_pt(&struc->cd.y);
  z = mod_float1_pt(&struc->cd.z);
  iresatm = mod_int1_pt(&struc->cd.iresatm);
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
        struct mdt_disulfide *bridge;
        g_array_set_size(ss, ss->len + 1);
        bridge = &g_array_index(ss, struct mdt_disulfide, ss->len - 1);
        bridge->atom1 = iat;
        bridge->atom2 = jat;
        bridge->res1 = iresatm[iat] - 1;
        bridge->res2 = iresatm[jat] - 1;
        mod_lognote("Disulfide bridge detected between residues %d and %d",
                    bridge->res1, bridge->res2);
      }
    }
  }
  disulfides = g_malloc(sizeof(struct mdt_disulfide_list));
  disulfides->nss = ss->len;
  disulfides->ss = (struct mdt_disulfide *)g_array_free(ss, FALSE);
  g_array_free(sg_atoms, TRUE);
  return disulfides;
}
