struct mod_mdt {
%immutable;
  int bin_type;
  int nelems;
  int nfeat;
};

struct mdt {
%immutable;
  struct mod_mdt base;
  gboolean symmetric;
  int n_proteins, n_protein_pairs;
  double sample_size;
  gboolean pdf;
};

struct mod_mdt_feature {
  int istart, iend, nbins, ifeat;
};

typedef enum {
  MOD_MDTB_FLOAT = 1,
  MOD_MDTB_DOUBLE,
  MOD_MDTB_INT32,
  MOD_MDTB_UINT32,
  MOD_MDTB_INT16,
  MOD_MDTB_UINT16,
  MOD_MDTB_INT8,
  MOD_MDTB_UINT8
} mod_mdt_bin_type;

struct mdt *mdt_new(mod_mdt_bin_type bin_type);
void mdt_free(struct mdt *mdt);

%inline %{
static struct mod_mdt_bin *mdt_library_bin_get(const struct mod_mdt *mdt,
                                               struct mdt_library *mlib,
                                               int nfeat, int nbin)
{
  int ifeat, istart;
  ifeat = mdt->features[nfeat].ifeat - 1;
  istart = mdt->features[nfeat].istart - 1;
  return &mlib->base.features[ifeat].bins[nbin + istart];
}

static struct mod_mdt_feature *mod_mdt_feature_get(const struct mod_mdt *mdt,
                                                   int nfeat)
{
  return &mdt->features[nfeat];
}
%}
