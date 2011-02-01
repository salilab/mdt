/** \file mdt_error.c      MDT error handling functions
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <glib.h>
#include "mdt_error.h"

/** Domain for MDT errors */
GQuark mdt_error_quark(void)
{
  return g_quark_from_static_string("mdt-error-quark");
}
