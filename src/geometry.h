/** \file geometry.h           Functions for calculating distances and angles.
 *
 *             Part of MDT, Copyright(c) 1989-2020 Andrej Sali
 */

#ifndef __MDT_GEOMETRY_H
#define __MDT_GEOMETRY_H

#include <glib.h>
#include "mdt_config.h"
#include "modeller.h"
#include "mdt_types.h"

G_BEGIN_DECLS

struct mdt_feature;

/** Return TRUE iff the coordinate value is undefined (-999.0) */
#define coordinate_undefined(x) (ABS(-999.0 - (x)) < 1e-4)

/** Return the squared distance between two coordinates.
    outrange is set to TRUE if the distance cannot be reliably calculated. */
MDTDLLLOCAL
float dist1sq(float x1, float y1, float z1, float x2, float y2, float z2,
              gboolean *outrange);

/** Return the bin index for the distance between two specified atoms in the
    same protein. */
MDTDLLLOCAL
int idist0(int ia1, int ia1p, const struct mod_structure *struc,
           const struct mdt_feature *feat);

/** Return the bin index for the angle between three specified atoms in the
    same protein. */
MDTDLLLOCAL
int iangle0(int ia1, int ia2, int ia3, const struct mod_structure *struc,
            const struct mdt_feature *feat);

/** Return the bin index for the dihedral angle between four specified atoms
    in the same protein. */
MDTDLLLOCAL
int idihedral0(int ia1, int ia2, int ia3, int ia4,
               const struct mod_structure *struc,
               const struct mdt_feature *feat);

/** Return the distance and the error on distance between two specified atoms
    in the same protein */
MDTDLLLOCAL
float dist0witherr(int ia1, int ia1p, const struct mod_structure *struc,
                   float *std, float errorscale);

/** Return the angle and the error on angle between three specified atoms
    in the same protein */
MDTDLLLOCAL
float angle0witherr(int ia1, int ia2, int ia3,
                    const struct mod_structure *struc, float *std,
                    float errorscale);

/** Return the dihedral and the error on dihedral between four specified atoms
    in the same protein */
MDTDLLLOCAL
float dihedral0witherr(int ia1, int ia2, int ia3, int ia4,
                       const struct mod_structure *struc, float *std,
                       float errorscale);

G_END_DECLS

#endif  /* __MDT_GEOMETRY_H */
