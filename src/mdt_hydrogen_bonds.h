/** \file mdt_hydrogen_bonds.h    Functions for handling hydrogen bonds
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_HBOND_H
#define __MDT_HBOND_H

#include <glib.h>
#include "mdt_config.h"
#include "mdt_types.h"
#include "mod_types.h"

G_BEGIN_DECLS

/** Return the indices of the "top-left" corner of the MDT. This must be freed
    by the user after use. */
MDTDLLLOCAL
int numb_hda(int ia, const int hb_iattyp[], const struct mod_coordinates *cd,
             const struct mdt_atom_class_list *atclass, float hbond_cutoff,
             int hbprop_type, int nbins);

/** Calculate H-bond protein satisfaction for the whole protein. */
MDTDLLLOCAL
float hb_satisfaction(const struct mod_coordinates *cd, const int iatta[],
                      const struct mdt_atom_class_list *atclass,
                      float hbond_cutoff);

G_END_DECLS

#endif  /* __MDT_HBOND_H */
