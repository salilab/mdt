/** \file mdt_geometry.h   Functions to calculate geometrical properties
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_GEOMETRY_H
#define __MDT_GEOMETRY_H

#include <glib.h>

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** Return the distance between two coordinates */
G_GNUC_INTERNAL
float dist1(float x1, float y1, float z1, float x2, float y2, float z2);

/** Return the angle between three coordinates */
G_GNUC_INTERNAL
float angle1(float x1, float y1, float z1, float x2, float y2, float z2,
             float x3, float y3, float z3);

/** Return the dihedral angle between four coordinates */
G_GNUC_INTERNAL
float dihedral1(float x1, float y1, float z1, float x2, float y2,
                float z2, float x3, float y3, float z3, float x4,
                float y4, float z4);

G_END_DECLS

#endif  /* __MDT_GEOMETRY_H */
