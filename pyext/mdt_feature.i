gboolean mdt_feature_periodic_get(const struct mdt_library *mlib, int ifeat);

void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins);

void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol);
