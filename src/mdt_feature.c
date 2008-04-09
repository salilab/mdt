/** \file mdt_feature.c    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"

/** Is the given feature type periodic? */
gboolean mdt_feature_is_periodic(int ifeat)
{
  switch (ifeat) {
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
  case 28:
  case 29:
  case 41:
  case 42:
  case 53:
  case 54:
  case 55:
  case 56:
  case 57:
  case 106:
  case 107:
  case 108:
  case 114:
    return TRUE;
  default:
    return FALSE;
  }
}
