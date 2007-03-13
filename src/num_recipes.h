/** \file num_recipes.h    Numerical functions, largely from Numerical Recipes.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __NUM_RECIPES_H
#define __NUM_RECIPES_H

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Return the number of combinations of n. */
G_GNUC_INTERNAL
int nperm(int n);

/** Return the complement of incomplete gamma function. */
G_GNUC_INTERNAL
float gammq(float a, float x);

#ifdef __cplusplus
}
#endif
#endif  /* __NUM_RECIPES_H */
