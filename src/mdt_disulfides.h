/** \file mdt_disulfides.h    Functions to handle disulfide bridges.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#ifndef __MDT_DISULFIDES_H
#define __MDT_DISULFIDES_H

#include <glib.h>
#include <stdio.h>
#include "mdt_config.h"
#include "modeller.h"
#include "mdt_types.h"

G_BEGIN_DECLS

struct mdt_properties;

/** Get list of disulfides. */
MDTDLLLOCAL
void get_disulfides(const struct mod_structure *struc,
                    const struct mod_sequence *seq,
                    const struct mod_libraries *libs,
                    struct mdt_properties *prop, int is);

G_END_DECLS

#endif  /* __MDT_DISULFIDES_H */
