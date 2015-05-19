/** \file num_recipes.h    Numerical functions, largely from Numerical Recipes.
 *
 *             Part of MDT, Copyright(c) 1989-2015 Andrej Sali
 */

#ifndef __NUM_RECIPES_H
#define __NUM_RECIPES_H

#include <glib.h>
#include "mdt_config.h"

G_BEGIN_DECLS

/** Return the number of combinations of n. */
MDTDLLLOCAL
int nperm(int n);

/** Return the complement of incomplete gamma function. */
MDTDLLLOCAL
float gammq(float a, float x);

G_END_DECLS

#endif  /* __NUM_RECIPES_H */
