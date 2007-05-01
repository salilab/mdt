struct mdt_bin {
%immutable;
  float rang1, rang2;
  char *symbol;
};

struct mod_mdt_library {
};

struct mdt_library {
%immutable;
  struct mod_mdt_library base;
%mutable;
  int deltai, deltaj;
  float hbond_cutoff;
  gboolean special_atoms, deltai_ali, deltaj_ali;
};

struct mdt_library *mdt_library_new(void);
void mdt_library_free(struct mdt_library *mlib);
gboolean mdt_atom_classes_read(const gchar *filename,
                               struct mdt_library *mlib, int natom,
                               GError **err);
gboolean mdt_hbond_read(const gchar *filename, struct mdt_library *mlib,
                        GError **err);
gboolean mdt_tuple_read(const gchar *filename, struct mdt_library *mlib, 
                        GError **err);
