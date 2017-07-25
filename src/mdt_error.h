/** \file mdt_error.h      MDT error handling.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
 */

#ifndef __MDT_ERROR_H
#define __MDT_ERROR_H

#include <glib.h>
#include "mdt_config.h"

G_BEGIN_DECLS

/** Domain for MDT errors */
#define MDT_ERROR mdt_error_quark()

/** MDT error types */
typedef enum {
  MDT_ERROR_IO,          /* Input/output error */
  MDT_ERROR_VALUE,       /* Incorrect value */
  MDT_ERROR_INDEX,       /* Index out of range */
  MDT_ERROR_FILE_FORMAT, /* File format error */
  MDT_ERROR_FAILED       /* Generic error */
} MDTError;

/** Domain for MDT errors */
MDTDLLEXPORT
GQuark mdt_error_quark(void);

G_END_DECLS

#endif  /* __MDT_ERROR_H */
