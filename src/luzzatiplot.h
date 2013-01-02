/** \file luzzatiplot.h       Estimate error from R-factor using Luzzati plots.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#ifndef __MDT_LUZZATI_PLOT_H
#define __MDT_LUZZATI_PLOT_H

#include <glib.h>
#include "mdt_config.h"

G_BEGIN_DECLS

/** Return error from R-factor using a pre-calculated Luzzati plot.
    The error returned is (error on coordinates)x(1/x-ray resolution)
 */
MDTDLLLOCAL
float luzzatiplot(int ind);

G_END_DECLS

#endif  /* __MDT_LUZZATI_PLOT_H */
