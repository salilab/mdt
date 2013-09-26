/** \file mdt_feature.h    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#ifndef __MDT_FEATURE_H
#define __MDT_FEATURE_H

#include <glib.h>
#include <mod_types.h>
#include "mdt_config.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** User-defined feature types */
typedef enum {
  MDT_FEATURE_NONE = 0,
  MDT_FEATURE_PROTEIN,
  MDT_FEATURE_PROTEIN_PAIR,
  MDT_FEATURE_RESIDUE,
  MDT_FEATURE_RESIDUE_PAIR,
  MDT_FEATURE_ALIGNED_RESIDUE,
  MDT_FEATURE_ALIGNED_RESIDUE_PAIR,
  MDT_FEATURE_ATOM,
  MDT_FEATURE_ATOM_PAIR,
  MDT_FEATURE_TUPLE,
  MDT_FEATURE_TUPLE_PAIR,
  MDT_FEATURE_BOND,
  MDT_FEATURE_GROUP
} mdt_feature_type;

struct mdt_properties;
struct mdt_library;
struct mdt_tuple;
struct mdt_bond;
struct mdt_feature;

typedef void (*mdt_cb_free)(void *data);

typedef int (*mdt_cb_feature_protein)(const struct mod_alignment *aln,
                                      int protein,
                                      struct mdt_properties *prop,
                                      const struct mdt_feature *feat,
                                      const struct mdt_library *mlib,
                                      const struct mod_libraries *libs,
                                      GError **err);

typedef int (*mdt_cb_feature_protein_pair)
    (const struct mod_alignment *aln, int protein1, int protein2,
     struct mdt_properties *prop,
     const struct mdt_feature *feat, const struct mdt_library *mlib,
     const struct mod_libraries *libs, GError **err);

typedef int (*mdt_cb_feature_residue)(const struct mod_alignment *aln,
                                      int protein, int residue,
                                      struct mdt_properties *prop,
                                      const struct mdt_feature *feat,
                                      const struct mdt_library *mlib,
                                      const struct mod_libraries *libs,
                                      GError **err);

typedef int (*mdt_cb_feature_residue_pair)
    (const struct mod_alignment *aln, int protein, int residue1, int residue2,
     struct mdt_properties *prop,
     const struct mdt_feature *feat, const struct mdt_library *mlib,
     const struct mod_libraries *libs, GError **err);

typedef int (*mdt_cb_feature_aligned_residue)
    (const struct mod_alignment *aln, int protein1, int protein2, int alnpos,
     struct mdt_properties *prop,
     const struct mdt_feature *feat, const struct mdt_library *mlib,
     const struct mod_libraries *libs, GError **err);

typedef int (*mdt_cb_feature_aligned_residue_pair)
    (const struct mod_alignment *aln, int protein1, int protein2, int alnpos1,
     int alnpos2, struct mdt_properties *prop,
     const struct mdt_feature *feat, const struct mdt_library *mlib,
     const struct mod_libraries *libs, GError **err);


typedef int (*mdt_cb_feature_atom)(const struct mod_alignment *aln,
                                   int protein, int atom,
                                   struct mdt_properties *prop,
                                   const struct mdt_feature *feat,
                                   const struct mdt_library *mlib,
                                   const struct mod_libraries *libs,
                                   GError **err);

typedef int (*mdt_cb_feature_atom_pair)(const struct mod_alignment *aln,
                                        int protein, int atom1, int atom2,
                                        struct mdt_properties *prop,
                                        const struct mdt_feature *feat,
                                        const struct mdt_library *mlib,
                                        const struct mod_libraries *libs,
                                        GError **err);

typedef int (*mdt_cb_feature_tuple)(const struct mod_alignment *aln,
                                    int protein, int atom,
                                    const struct mdt_tuple *tuple,
                                    struct mdt_properties *prop,
                                    const struct mdt_feature *feat,
                                    const struct mdt_library *mlib,
                                    const struct mod_libraries *libs,
                                    GError **err);

typedef int (*mdt_cb_feature_tuple_pair)(const struct mod_alignment *aln,
                                         int protein, int atom1,
                                         const struct mdt_tuple *tuple1,
                                         int atom2,
                                         const struct mdt_tuple *tuple2,
                                         struct mdt_properties *prop,
                                         const struct mdt_feature *feat,
                                         const struct mdt_library *mlib,
                                         const struct mod_libraries *libs,
                                         GError **err);

typedef int (*mdt_cb_feature_bond)(const struct mod_alignment *aln,
                                   int protein, const struct mdt_bond *bond,
                                   struct mdt_properties *prop,
                                   const struct mdt_feature *feat,
                                   const struct mdt_library *mlib,
                                   const struct mod_libraries *libs,
                                   GError **err);

typedef int (*mdt_cb_feature_group)(int bin1, int bin2,
                                   struct mdt_properties *prop,
                                   const struct mdt_feature *feat,
                                   const struct mdt_library *mlib,
                                   const struct mod_libraries *libs,
                                   GError **err);

/** User-defined protein feature */
struct mdt_feature_protein {
  int protein;
  mdt_cb_feature_protein getbin;
};

/** User-defined protein pair feature */
struct mdt_feature_protein_pair {
  int protein1;
  int protein2;
  mdt_cb_feature_protein_pair getbin;
};

/** User-defined residue feature */
struct mdt_feature_residue {
  int protein;
  int delta;
  int align_delta;
  gboolean pos2;
  int bin_seq_outrange;
  mdt_cb_feature_residue getbin;
};

/** User-defined residue pair feature */
struct mdt_feature_residue_pair {
  int protein;
  mdt_cb_feature_residue_pair getbin;
};

/** User-defined aligned residue feature */
struct mdt_feature_aligned_residue {
  int protein1;
  int protein2;
  mdt_cb_feature_aligned_residue getbin;
};

/** User-defined aligned residue pair feature */
struct mdt_feature_aligned_residue_pair {
  int protein1;
  int protein2;
  mdt_cb_feature_aligned_residue_pair getbin;
};

/** User-defined atom feature */
struct mdt_feature_atom {
  gboolean pos2;
  mdt_cb_feature_atom getbin;
};

/** User-defined atom pair feature */
struct mdt_feature_atom_pair {
  mdt_cb_feature_atom_pair getbin;
};

/** User-defined tuple feature */
struct mdt_feature_tuple {
  gboolean pos2;
  mdt_cb_feature_tuple getbin;
};

/** User-defined tuple pair feature */
struct mdt_feature_tuple_pair {
  mdt_cb_feature_tuple_pair getbin;
};

/** User-defined tuple pair feature */
struct mdt_feature_bond {
  int type;
  mdt_cb_feature_bond getbin;
};

/** User-defined group feature */
struct mdt_feature_group {
  int ifeat1, ifeat2;
  mdt_cb_feature_group getbin;
};

/** MDT feature */
struct mdt_feature {
  /** Base Modeller feature (stores bins, name, Modeller properties such as
      accessible surface area to precalculate). Only valid during a scan. */
  struct mod_mdt_libfeature *base;
  mdt_feature_type type;
  union {
    struct mdt_feature_protein protein;
    struct mdt_feature_protein_pair protein_pair;
    struct mdt_feature_residue residue;
    struct mdt_feature_residue_pair residue_pair;
    struct mdt_feature_aligned_residue_pair aligned_residue_pair;
    struct mdt_feature_aligned_residue aligned_residue;
    struct mdt_feature_atom atom;
    struct mdt_feature_atom_pair atom_pair;
    struct mdt_feature_tuple tuple;
    struct mdt_feature_tuple_pair tuple_pair;
    struct mdt_feature_bond bond;
    struct mdt_feature_group group;
  } u;
  void *data;
  mdt_cb_free freefunc;
  /** TRUE if the feature range is periodic (e.g. for a dihedral) */
  gboolean periodic;
  /** TRUE during scans if the bins are of uniform width */
  gboolean uniform_bins;
  /** If uniform bins, the reciprocal of the bin width */
  float inverse_bin_width;
  /** The square of the largest value that can be binned during a scan */
  float max_range_squared;
};

/** Get the index of the undefined bin. */
#define mdt_feature_undefined_bin_get(feat) ((feat)->base->nbins)

/** Add a protein feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_protein_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein,
                            mdt_cb_feature_protein getbin, void *data,
                            mdt_cb_free freefunc, GError **err);

/** Add a protein pair feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_protein_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein1,
                                 int protein2,
                                 mdt_cb_feature_protein_pair getbin, void *data,
                                 mdt_cb_free freefunc, GError **err);

/** Add a residue feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_residue_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein, int delta,
                            int align_delta, gboolean pos2,
                            int bin_seq_outrange, mdt_cb_feature_residue getbin,
                            void *data, mdt_cb_free freefunc, GError **err);

/** Add a residue pair feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_residue_pair_add(struct mdt_library *mlib, const char *name,
                                 mod_mdt_calc precalc_type, int protein,
                                 gboolean asymmetric,
                                 mdt_cb_feature_residue_pair getbin, void *data,
                                 mdt_cb_free freefunc, GError **err);

/** Add an aligned residue feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_aligned_residue_add(struct mdt_library *mlib, const char *name,
                                    mod_mdt_calc precalc_type, int protein1,
                                    int protein2,
                                    mdt_cb_feature_aligned_residue getbin,
                                    void *data, mdt_cb_free freefunc,
                                    GError **err);

/** Add an aligned residue pair feature.
    \return the index of the new feature, or -1 on error. */
MDTDLLEXPORT
int mdt_feature_aligned_residue_pair_add(
    struct mdt_library *mlib, const char *name, mod_mdt_calc precalc_type,
    int protein1, int protein2, gboolean asymmetric,
    mdt_cb_feature_aligned_residue_pair getbin, void *data,
    mdt_cb_free freefunc, GError **err);

/** Add an atom feature.
    \note The system is automatically instructed to read in PDB files.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_atom_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type, gboolean pos2,
                         mdt_cb_feature_atom getbin, void *data,
                         mdt_cb_free freefunc);

/** Add an atom pair feature.
    \note The system is automatically instructed to read in PDB files.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_atom_pair_add(struct mdt_library *mlib, const char *name,
                              mod_mdt_calc precalc_type, gboolean asymmetric,
                              mdt_cb_feature_atom_pair getbin, void *data,
                              mdt_cb_free freefunc);

/** Add a tuple feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_tuple_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, gboolean pos2,
                          mdt_cb_feature_tuple getbin, void *data,
                          mdt_cb_free freefunc);

/** Add a tuple pair feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_tuple_pair_add(struct mdt_library *mlib, const char *name,
                               mod_mdt_calc precalc_type,
                               mdt_cb_feature_tuple_pair getbin, void *data,
                               mdt_cb_free freefunc);

/** Add a bond feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_bond_add(struct mdt_library *mlib, const char *name,
                         mod_mdt_calc precalc_type,
                         mdt_cb_feature_bond getbin, void *data,
                         mdt_cb_free freefunc);

/** Add an angle feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_angle_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type,
                          mdt_cb_feature_bond getbin, void *data,
                          mdt_cb_free freefunc);

/** Add a dihedral feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_dihedral_add(struct mdt_library *mlib, const char *name,
                             mod_mdt_calc precalc_type,
                             mdt_cb_feature_bond getbin, void *data,
                             mdt_cb_free freefunc);

/** Add a group feature.
    \return the index of the new feature */
MDTDLLEXPORT
int mdt_feature_group_add(struct mdt_library *mlib, const char *name,
                          mod_mdt_calc precalc_type, int ifeat1, int ifeat2,
                          mdt_cb_feature_group getbin, void *data,
                          mdt_cb_free freefunc);

/** Add a data file type needed for a given feature */
MDTDLLEXPORT
void mdt_feature_add_needed_file(struct mdt_library *mlib, int ifeat,
                                 mod_mdt_file filetype);

/** Set whether a feature's range is periodic (FALSE by default) */
MDTDLLEXPORT
void mdt_feature_periodic_set(struct mdt_library *mlib, int ifeat,
                              gboolean periodic);

/** Get whether a feature's range is periodic */
MDTDLLEXPORT
gboolean mdt_feature_periodic_get(const struct mdt_library *mlib, int ifeat);

/** Set the number of bins for a given feature.
    \note One extra bin is always created - the 'undefined' bin. */
MDTDLLEXPORT
void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins);

/** Set the range and symbol for a given bin and feature. */
MDTDLLEXPORT
void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol);

G_END_DECLS

#endif  /* __MDT_FEATURE_H */
