struct mdt_type {
%immutable;
  int nelems;
  int nfeat;
  int n_proteins, n_protein_pairs;
  double sample_size;
  gboolean pdf;
  struct f_int1_array istart, iend, nbins, ifeat;
  struct f_float1_array bin;
};

%inline %{
static struct mdt_bin *mdt_library_bin_get(const struct mdt_type *mdt,
                                           struct mdt_library *mlib, int nfeat,
                                           int nbin)
{
  int ifeat, istart;
  ifeat = f_int1_get(&mdt->ifeat, nfeat) - 1;
  istart = f_int1_get(&mdt->istart, nfeat) - 1;
  return &mlib->base.features[ifeat].bins[nbin + istart];
}
%}
