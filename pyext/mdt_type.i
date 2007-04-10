struct mdt_type {
%immutable;
  int nelems;
  int nfeat;
  int n_proteins, n_protein_pairs;
  double sample_size;
  gboolean pdf;
  struct f_float1_array bin;
};

struct mdt_feature {
  int istart, iend, nbins, ifeat;
};

%inline %{
static struct mdt_bin *mdt_library_bin_get(const struct mdt_type *mdt,
                                           struct mdt_library *mlib, int nfeat,
                                           int nbin)
{
  int ifeat, istart;
  ifeat = mdt->features[nfeat].ifeat - 1;
  istart = mdt->features[nfeat].istart - 1;
  return &mlib->base.features[ifeat].bins[nbin + istart];
}

static struct mdt_feature *mdt_type_feature_get(const struct mdt_type *mdt,
                                                int nfeat)
{
  return &mdt->features[nfeat];
}
%}
