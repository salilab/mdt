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
  gboolean deltai_ali, deltaj_ali;
};

struct mdt_library *mdt_library_new(void);
void mdt_library_free(struct mdt_library *mlib);
