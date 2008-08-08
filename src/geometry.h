/** \file geometry.h           Functions for calculating distances and angles.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_GEOMETRY_H
#define __MDT_GEOMETRY_H

#include <glib.h>
#include "mdt_config.h"
#include "modeller.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** Return the bin index for the distance between two specified atoms in the
    same protein. */
MDTDLLLOCAL
int idist0(int ia1, int ia1p, const struct mod_structure *struc,
           const struct mod_mdt_libfeature *feat);

/** Return the bin index for the angle between three specified atoms in the
    same protein. */
MDTDLLLOCAL
int iangle0(int ia1, int ia2, int ia3, const struct mod_structure *struc,
            const struct mod_mdt_libfeature *feat);

/** Return the bin index for the dihedral angle between four specified atoms
    in the same protein. */
MDTDLLLOCAL
int idihedral0(int ia1, int ia2, int ia3, int ia4,
               const struct mod_structure *struc,
               const struct mod_mdt_libfeature *feat);

G_END_DECLS

#endif  /* __MDT_GEOMETRY_H */
