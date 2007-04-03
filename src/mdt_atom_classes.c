/** \file mdt_atom_classes.c   Functions to handle atom classes.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdio.h>
#include <glib.h>
#include "mod_file.h"
#include "util.h"
#include "mdt_error.h"
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
  for (i = 0; i < natom; i++) {
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
                mymsg);
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

static gboolean read_atmgrp(GScanner *scanner, GArray *classes, GArray **types,
                            gboolean read_hbond, GError **err)
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
                          GArray *classes, GArray *types, GError **err)
{
  int i;
  struct mdt_atom_type atype;
  gboolean retval = TRUE;
  if (classes->len == 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "ATOM before ATMGRP");
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
                                       gboolean read_hbond, GError **err)
{
  static const char *grpnames[] = { "ATMGRP", "BNDGRP", "ANGGRP", "DIHGRP" };
  static const char *atnames[] = { "ATOM", "BOND", "ANGLE", "DIHEDRAL" };
  const char *sym[2];
  gboolean retval = TRUE;
  GScanner *scanner = g_scanner_new(NULL);
  GArray *classes = g_array_new(FALSE, FALSE, sizeof(struct mdt_atom_class));
  GArray *types = g_array_new(FALSE, FALSE, sizeof(struct mdt_atom_type));

  sym[0] = grpnames[atclass->natom - 1];
  sym[1] = atnames[atclass->natom - 1];
  g_scanner_add_symbol(scanner, sym[0], GINT_TO_POINTER(0));
  g_scanner_add_symbol(scanner, sym[1], GINT_TO_POINTER(1));
  scanner->input_name = filename;
  scanner->config->int_2_float = TRUE;
  g_scanner_input_text(scanner, text, filelen);
  while (g_scanner_get_next_token(scanner) != G_TOKEN_EOF && retval) {
    if (scanner->token == G_TOKEN_SYMBOL) {
      if (GPOINTER_TO_INT(scanner->value.v_symbol) == 0) {
        retval = read_atmgrp(scanner, classes, &types, read_hbond, err);
      } else {
        retval = read_atom(scanner, atclass, classes, types, err);
      }
    } else {
      mod_g_scanner_unexp(scanner, G_TOKEN_SYMBOL, NULL, "ATOM or ATMGRP",
                          err);
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

static void update_mdt_feat_atclass(struct mdt_feature *feat,
                                    const struct mdt_atom_class_list *atclass,
                                    int natom)
{
  int i;
  mdt_feature_nbins_set(feat, atclass->nclass + 1);
  for (i = 0; i < atclass->nclass; i++) {
    g_free(feat->bins[i].symbol);
    feat->bins[i].symbol = g_strdup(atclass->classes[i].name);
  }
  g_free(feat->bins[atclass->nclass].symbol);
  feat->bins[atclass->nclass].symbol = g_strdup("U");
}

static gboolean read_atom_class_file(const gchar *filename,
                                     struct mdt_library *mlib,
                                     struct mdt_atom_class_list *atclass,
                                     gboolean read_hbond, GError **err)
{
  struct mod_file file_info;
  FILE *fp;
  char *text;
  int ierr;
  gboolean retval;
  unsigned filelen;

  fp = mdt_open_file(filename, "r", &file_info, err);
  if (!fp) {
    return FALSE;
  }

  mod_file_read_contents(fp, &text, &filelen, &ierr);
  if (ierr != 0) {
    handle_modeller_error(err);
    return FALSE;
  } else {
    retval = scan_atom_classes_file(filename, text, filelen, atclass,
                                    read_hbond, err);
    g_free(text);
  }

  return mdt_close_file(fp, &file_info, err) && retval;
}

/** Read atom class information from a file; return TRUE on success. */
gboolean mdt_atom_classes_read(const gchar *filename,
                               struct mdt_library *mlib, int natom,
                               GError **err)
{
  gboolean retval;
  retval = read_atom_class_file(filename, mlib, mlib->atclass[natom - 1],
                                FALSE, err);
  if (retval) {
    int ifeat;
    struct mdt_feature *feat;

    /* MDT features; 79 = atom, 109 = bond, 111 = angle, 113 = dihedral types */
    if (natom == 1) {
      ifeat = 78;
    } else {
      ifeat = 104 + natom * 2;
    }

    /* Set MDT symbols */
    feat = &mlib->base.features[ifeat];
    update_mdt_feat_atclass(feat, mlib->atclass[natom - 1], natom);
  }
  return retval;
}

/** Read hydrogen bond class information from a file; return TRUE on success. */
gboolean mdt_hbond_read(const gchar *filename, struct mdt_library *mlib, 
                        GError **err)
{
  return read_atom_class_file(filename, mlib, mlib->hbond, TRUE, err);
}
