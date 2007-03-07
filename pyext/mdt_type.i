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
