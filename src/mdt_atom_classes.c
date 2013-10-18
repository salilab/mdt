/** \file mdt_atom_classes.c   Functions to handle atom classes.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <stdio.h>
#include <glib.h>
#include "mod_file.h"
#include "util.h"
#include "mdt_error.h"
#include "mdt_hdf5.h"
#include "mdt_atom_classes.h"

/** Make a new atom class list. */
struct mdt_atom_class_list *mdt_atom_class_list_new(int natom)
{
  struct mdt_atom_class_list *c = g_malloc(sizeof(struct mdt_atom_class_list));
  c->natom = natom;
  c->nclass = 0;
  c->classes = NULL;
  return c;
}

static void mdt_atom_type_free(struct mdt_atom_type *attype, int natom)
{
  int i;
  for (i = 0; i <= natom; i++) {
    g_free(attype->names[i]);
  }
  g_free(attype->names);
}

static void mdt_atom_class_free(struct mdt_atom_class *atclass, int natom)
{
  int i;
  for (i = 0; i < atclass->ntypes; i++) {
    mdt_atom_type_free(&atclass->types[i], natom);
  }
  g_free(atclass->types);
  g_free(atclass->name);
}

/** Free an existing atom class list. */
void mdt_atom_class_list_free(struct mdt_atom_class_list *atclass)
{
  int i;
  for (i = 0; i < atclass->nclass; i++) {
    mdt_atom_class_free(&atclass->classes[i], atclass->natom);
  }
  g_free(atclass->classes);
}

static void my_handler(GScanner *scanner, gchar *message, gboolean error)
{
  char *mymsg;
  mymsg = g_strdup_printf("%s:%d:%d: %s", scanner->input_name,
                          g_scanner_cur_line(scanner),
                          g_scanner_cur_position(scanner), message);
  puts(mymsg);
  if (error) {
    g_set_error((GError **)scanner->user_data, MDT_ERROR, MDT_ERROR_FAILED,
                "%s", mymsg);
  }
  g_free(mymsg);
}

static void mod_g_scanner_unexp(GScanner *scanner, GTokenType expected_token,
                                const gchar *identifier_spec,
                                const gchar *symbol_spec, GError **err)
{
  scanner->msg_handler = my_handler;
  scanner->user_data = err;
  g_scanner_unexp_token(scanner, expected_token, identifier_spec, symbol_spec,
                        NULL, NULL, TRUE);
}

static gboolean read_atmgrp_hbond(GScanner *scanner,
                                  struct mdt_atom_class *aclass, GError **err)
{
  gboolean retval = TRUE;
  int i;
  for (i = 0; i < 3 && retval; i++) {
    float sign;
    GTokenType token = g_scanner_get_next_token(scanner);
    if (token == '-' || token == '+') {
      sign = (token == '+' ? 1.0 : -1.0);
      token = g_scanner_get_next_token(scanner);
    } else {
      sign = 1.0;
    }
    if (token == G_TOKEN_FLOAT) {
      aclass->hb_property[i] = sign * scanner->value.v_float;
    } else {
      mod_g_scanner_unexp(scanner, G_TOKEN_FLOAT, NULL, NULL, err);
      retval = FALSE;
    }
  }
  return retval;
}

static gboolean read_atmgrp(GScanner *scanner, GArray *classes,
                            GArray **types, gboolean read_hbond, GError **err)
{
  if (g_scanner_get_next_token(scanner) == G_TOKEN_STRING) {
    struct mdt_atom_class newclass;
    if (classes->len > 0) {
      struct mdt_atom_class *c = &g_array_index(classes, struct mdt_atom_class,
                                                classes->len - 1);
      c->ntypes = (*types)->len;
      c->types = (struct mdt_atom_type *)g_array_free(*types, FALSE);
      *types = g_array_new(FALSE, FALSE, sizeof(struct mdt_atom_type));
    }
    newclass.name = g_strdup(scanner->value.v_string);
    if (read_hbond && !read_atmgrp_hbond(scanner, &newclass, err)) {
      g_free(newclass.name);
      return FALSE;
    }
    g_array_append_val(classes, newclass);
    return TRUE;
  } else {
    mod_g_scanner_unexp(scanner, G_TOKEN_STRING, NULL, NULL, err);
    return FALSE;
  }
}

static gboolean read_atom(GScanner *scanner,
                          struct mdt_atom_class_list *atclass,
                          GArray *classes, GArray *types, const char *sym[2],
                          GError **err)
{
  int i;
  struct mdt_atom_type atype;
  gboolean retval = TRUE;
  if (classes->len == 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "%s:%d:%d: %s line encountered without a preceding %s line",
                scanner->input_name, g_scanner_cur_line(scanner),
                g_scanner_cur_position(scanner), sym[1], sym[0]);
    return FALSE;
  }

  atype.names = g_malloc(sizeof(char *) * (atclass->natom + 1));
  for (i = 0; i < atclass->natom + 1 && retval; i++) {
    if (g_scanner_get_next_token(scanner) == G_TOKEN_STRING) {
      atype.names[i] = g_strdup(scanner->value.v_string);
    } else {
      int j;
      for (j = 0; j < i; j++) {
        g_free(atype.names[j]);
      }
      mod_g_scanner_unexp(scanner, G_TOKEN_STRING, NULL, NULL, err);
      retval = FALSE;
    }
  }

  if (retval) {
    g_array_append_val(types, atype);
  } else {
    g_free(atype.names);
  }
  return retval;
}

static gboolean scan_atom_classes_file(const char *filename, const char *text,
                                       unsigned filelen,
                                       struct mdt_atom_class_list *atclass,
                                       gboolean read_hbond, gboolean tuples,
                                       GError **err)
{
  static const char *grpnames[] = { "ATMGRP", "BNDGRP", "ANGGRP", "DIHGRP" };
  static const char *atnames[] = { "ATOM", "BOND", "ANGLE", "DIHEDRAL" };
  const char *sym[2];
  gboolean retval = TRUE;
  GScanner *scanner = g_scanner_new(NULL);
  GArray *classes = g_array_new(FALSE, FALSE, sizeof(struct mdt_atom_class));
  GArray *types = g_array_new(FALSE, FALSE, sizeof(struct mdt_atom_type));

  if (tuples) {
    sym[0] = "DBLGRP";
    sym[1] = "TRPGRP";
    g_scanner_add_symbol(scanner, sym[0], GINT_TO_POINTER(2));
    g_scanner_add_symbol(scanner, sym[1], GINT_TO_POINTER(3));
  } else {
    sym[0] = grpnames[atclass->natom - 1];
    sym[1] = atnames[atclass->natom - 1];
    g_scanner_add_symbol(scanner, sym[0], GINT_TO_POINTER(0));
    g_scanner_add_symbol(scanner, sym[1], GINT_TO_POINTER(1));
  }
  scanner->input_name = filename;
  scanner->config->int_2_float = TRUE;
  g_scanner_input_text(scanner, text, filelen);
  while (g_scanner_get_next_token(scanner) != G_TOKEN_EOF && retval) {
    if (scanner->token == G_TOKEN_SYMBOL) {
      int symb = GPOINTER_TO_INT(scanner->value.v_symbol);
      switch (symb) {
      case 0:
        retval = read_atmgrp(scanner, classes, &types, read_hbond, err);
        break;
      case 1:
        retval = read_atom(scanner, atclass, classes, types, sym, err);
        break;
      case 2:
      case 3:
        atclass->natom = symb;
        g_scanner_remove_symbol(scanner, sym[0]);
        g_scanner_remove_symbol(scanner, sym[1]);
        if (symb == 2) {
          sym[0] = "DBLGRP";
          sym[1] = "DOUBLET";
        } else {
          sym[0] = "TRPGRP";
          sym[1] = "TRIPLET";
        }
        g_scanner_add_symbol(scanner, sym[0], GINT_TO_POINTER(0));
        g_scanner_add_symbol(scanner, sym[1], GINT_TO_POINTER(1));
        retval = read_atmgrp(scanner, classes, &types, read_hbond, err);
      }
    } else {
      char *symbols = g_strdup_printf("%s or %s", sym[0], sym[1]);
      mod_g_scanner_unexp(scanner, G_TOKEN_SYMBOL, NULL, symbols, err);
      g_free(symbols);
      retval = FALSE;
    }
  }
  atclass->nclass = classes->len;
  atclass->classes = (struct mdt_atom_class *)g_array_free(classes, FALSE);
  if (atclass->nclass > 0) {
    struct mdt_atom_class *c = &atclass->classes[atclass->nclass - 1];
    c->ntypes = types->len;
    c->types = (struct mdt_atom_type *)g_array_free(types, FALSE);
  }
  g_scanner_destroy(scanner);
  return retval;
}

/** Set the number of bins and the bin symbols for atom class features */
void update_mdt_feat_atclass(struct mod_mdt_libfeature *feat,
                             const struct mdt_atom_class_list *atclass)
{
  int i;
  mod_mdt_libfeature_nbins_set(feat, atclass->nclass + 1);
  for (i = 0; i < atclass->nclass; i++) {
    g_free(feat->bins[i].symbol);
    feat->bins[i].symbol = g_strdup(atclass->classes[i].name);
    feat->bins[i].rang1 = i;
    feat->bins[i].rang2 = i + 1;
  }
  g_free(feat->bins[atclass->nclass].symbol);
  feat->bins[atclass->nclass].symbol = g_strdup("U");
}

static gboolean write_atom_class_info(hid_t loc_id,
                                      const struct mdt_atom_class_list *atclass,
                                      gboolean write_hbond)
{
  int i, j, k;
  hsize_t featdim[2], typedim[2];
  hid_t dtype;
  gboolean retval;
  int *ntypes, total_ntypes;
  float *hbond = NULL;
  char **class_name, **type_name, **pt;

  ntypes = g_malloc(atclass->nclass * sizeof(int));
  class_name = g_malloc(atclass->nclass * sizeof(char *));

  if (write_hbond) {
    float *pt = hbond = g_malloc(atclass->nclass * 3 * sizeof(float));
    for (i = 0; i < atclass->nclass; ++i) {
      struct mdt_atom_class *cls = &atclass->classes[i];
      *(pt++) = cls->hb_property[0];
      *(pt++) = cls->hb_property[1];
      *(pt++) = cls->hb_property[2];
    }
  }

  total_ntypes = 0;
  for (i = 0; i < atclass->nclass; ++i) {
    ntypes[i] = atclass->classes[i].ntypes;
    total_ntypes += ntypes[i];
    class_name[i] = atclass->classes[i].name;
  }
  type_name = g_malloc(total_ntypes * (atclass->natom + 1) * sizeof(char *));
  typedim[0] = total_ntypes;
  typedim[1] = atclass->natom + 1;
  pt = type_name;
  for (i = 0; i < atclass->nclass; ++i) {
    struct mdt_atom_class *cls = &atclass->classes[i];
    for (j = 0; j < cls->ntypes; ++j) {
      struct mdt_atom_type *typ = &cls->types[j];
      for (k = 0; k < typedim[1]; ++k) {
        *(pt++) = typ->names[k];
      }
    }
  }

  featdim[0] = atclass->nclass;
  retval = (dtype = H5Tcopy(H5T_C_S1)) >= 0
           && H5Tset_size(dtype, H5T_VARIABLE) >= 0
           && H5Tset_strpad(dtype, H5T_STR_NULLTERM) >= 0
           && H5LTmake_dataset(loc_id, "class_names", 1, featdim, dtype,
                               class_name) >= 0
           && H5LTmake_dataset_int(loc_id, "ntypes", 1, featdim, ntypes) >= 0
           && H5LTmake_dataset(loc_id, "type_names", 2, typedim, dtype,
                               type_name) >= 0
           && H5Tclose(dtype) >= 0;

  if (write_hbond) {
    featdim[1] = 3;
    retval = retval && H5LTmake_dataset_float(loc_id, "hb_property", 2,
                                              featdim, hbond) >= 0;
    g_free(hbond);
  }

  g_free(ntypes);
  g_free(class_name);
  g_free(type_name);
  return retval;
}

static gboolean write_atom_class_file(hid_t loc_id, const char *name,
                                      const struct mdt_library *mlib,
                                      const struct mdt_atom_class_list *atclass,
                                      gboolean write_hbond)
{
  hid_t group_id;

  group_id = H5Gcreate(loc_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (group_id < 0) {
    return FALSE;
  }

  return write_atom_class_info(group_id, atclass, write_hbond)
         && H5Gclose(group_id) >= 0;
}

static gboolean read_atom_class_file(const gchar *filename,
                                     struct mdt_library *mlib,
                                     struct mdt_atom_class_list *atclass,
                                     gboolean read_hbond, gboolean tuples,
                                     GError **err)
{
  struct mod_file *fh;
  char *text;
  gboolean retval;
  unsigned filelen;

  if (mlib->feature_added) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "You must read all atom classes BEFORE creating any features");
    return FALSE;
  }
  fh = mdt_open_file(filename, "r", err);
  if (!fh) {
    return FALSE;
  }

  if (!mod_file_read_contents(fh, &text, &filelen)) {
    handle_modeller_error(err);
    retval = FALSE;
  } else {
    retval = scan_atom_classes_file(filename, text, filelen, atclass,
                                    read_hbond, tuples, err);
    g_free(text);
  }

  return mdt_close_file(fh, err) && retval;
}

/** Read atom class information from a file; return TRUE on success. */
gboolean mdt_atom_classes_read(const gchar *filename,
                               struct mdt_library *mlib, int natom,
                               GError **err)
{
  return read_atom_class_file(filename, mlib, mlib->atclass[natom - 1],
                              FALSE, FALSE, err);
}

/** Read hydrogen bond class information from a file; return TRUE on success. */
gboolean mdt_hbond_read(const gchar *filename, struct mdt_library *mlib,
                        GError **err)
{
  return read_atom_class_file(filename, mlib, mlib->hbond, TRUE, FALSE, err);
}

/** Read tuple class information from a file; return TRUE on success. */
gboolean mdt_tuple_read(const gchar *filename, struct mdt_library *mlib,
                        GError **err)
{
  return read_atom_class_file(filename, mlib, mlib->tupclass, FALSE, TRUE, err);
}

gboolean mdt_tuple_write(hid_t loc_id, const struct mdt_library *mlib)
{
  return write_atom_class_file(loc_id, "tuples", mlib, mlib->tupclass,
                               FALSE);
}

gboolean mdt_atom_class_write(hid_t loc_id, const struct mdt_library *mlib)
{
  return write_atom_class_file(loc_id, "atom_classes", mlib, mlib->atclass[0],
                               FALSE);
}

gboolean mdt_bond_class_write(hid_t loc_id, const struct mdt_library *mlib)
{
  return write_atom_class_file(loc_id, "bond_classes", mlib, mlib->atclass[1],
                               FALSE);
}

gboolean mdt_angle_class_write(hid_t loc_id, const struct mdt_library *mlib)
{
  return write_atom_class_file(loc_id, "angle_classes", mlib, mlib->atclass[2],
                               FALSE);
}

gboolean mdt_dihedral_class_write(hid_t loc_id, const struct mdt_library *mlib)
{
  return write_atom_class_file(loc_id, "dihedral_classes", mlib,
                               mlib->atclass[3], FALSE);
}

gboolean mdt_hbond_write(hid_t loc_id, const struct mdt_library *mlib)
{
  return write_atom_class_file(loc_id, "hbond", mlib, mlib->hbond, TRUE);
}
