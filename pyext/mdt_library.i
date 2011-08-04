struct mod_mdt_bin {
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
  float hbond_cutoff;
  gboolean special_atoms;
};

struct mdt_library *mdt_library_new(const struct mod_libraries *libs);
void mdt_library_free(struct mdt_library *mlib);
gboolean mdt_atom_classes_read(const gchar *filename,
                               struct mdt_library *mlib, int natom,
                               GError **err);
gboolean mdt_hbond_read(const gchar *filename, struct mdt_library *mlib,
                        GError **err);
gboolean mdt_tuple_read(const gchar *filename, struct mdt_library *mlib, 
                        GError **err);

%inline %{
static void mdt_library_distance_atoms_set(struct mdt_library *mlib,
                                           char *dstatm1, char *dstatm2)
{
  g_free(mlib->distance_atoms[0]);
  g_free(mlib->distance_atoms[1]);
  mlib->distance_atoms[0] = g_strdup(dstatm1);
  mlib->distance_atoms[1] = g_strdup(dstatm2);
}
%}
