/** \file mdt.h            Functions to handle MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#ifndef __MDT_MDT_H
#define __MDT_MDT_H

#include <glib.h>
#include "mdt_config.h"
#include "mdt_error.h"
#include "mdt_types.h"
#include "mdt_atom_classes.h"
#include "mod_types.h"

#include "mdt_alignment.h"

G_BEGIN_DECLS

/** Smooth a histogram or the 2D plot with a uniform prior.
    Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_smooth(const struct mdt *mdtin, struct mdt *mdtout,
                    int dimensions, float weight, GError **err);

/** Transform an MDT with an exp function */
MDTDLLEXPORT
void mdt_exp_transform(struct mod_mdt *mdt, float offset, float expoffset,
                       float multiplier, float power);

/** Transform an MDT with a linear function */
MDTDLLEXPORT
void mdt_linear_transform(struct mod_mdt *mdt, float offset, float multiplier);

/** Transform an MDT with an inverse function */
MDTDLLEXPORT
void mdt_inverse_transform(struct mod_mdt *mdt, float offset, float multiplier,
                           float undefined);

/** Transform an MDT with a log function */
MDTDLLEXPORT
void mdt_log_transform(struct mod_mdt *mdt, float offset, float multiplier,
                       float undefined);

/** Offset an MDT by the minimum value. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_offset_min(struct mod_mdt *mdt, int dimensions, GError **err);

/** Close an MDT so that it is useful for creating periodic splines.
    Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_close(struct mod_mdt *mdt, int dimensions, GError **err);

/** Get the entropy of the dependent variable. */
MDTDLLEXPORT
float mdt_entropy_hx(const struct mod_mdt *mdt, GError **err);

/** Write chi^2, entropies and dependencies of pdf p(x/y,z,...) to the log.
    Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_entropy_full(const struct mod_mdt *mdt,
                          const struct mdt_library *mlib, GError **err);

/** Read in an MDT in text format. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_read(struct mdt *mdt, const struct mdt_library *mlib,
                  const char *filename, GError **err);

/** Write out an MDT in text format. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_write(const struct mdt *mdt, const struct mdt_library *mlib,
                   const char *filename, gboolean write_preamble, GError **err);

/** Write out an MDT in HDF5 format. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_write_hdf5(const struct mdt *mdt, const struct mdt_library *mlib,
                        const char *filename, int gzip, int chunk_size,
                        GError **err);

/** Read in an MDT in HDF5 format. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_read_hdf5(struct mdt *mdt, const struct mdt_library *mlib,
                       const char *filename, GError **err);

/** Normalize an MDT. */
MDTDLLEXPORT
gboolean mdt_normalize(const struct mdt *mdtin, struct mdt *mdtout,
                       const struct mdt_library *mlib, int dimensions,
                       const float dx_dy[], int n_dx_dy, gboolean to_zero,
                       gboolean to_pdf, GError **err);

/** Integrate an MDT. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_integrate(const struct mdt *mdtin, struct mdt *mdtout,
                       const int features[], int n_features, GError **err);

/** Reshape an MDT. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_reshape(const struct mdt *mdtin, struct mdt *mdtout,
                     const int features[], int n_features, const int offset[],
                     int n_offset, const int shape[], int n_shape,
                     GError **err);

/** Clear an MDT (set all elements to zero). */
MDTDLLEXPORT
void mdt_clear(struct mdt *mdt);

/** Clear the MDT array, and set feature types. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_make(struct mdt *mdt, const struct mdt_library *mlib,
                  const int features[], int n_features, const int shape[],
                  int n_shape, GError **err);

/** Make mdtout a copy of mdtin. */
MDTDLLEXPORT
void mdt_copy(const struct mdt *mdtin, struct mdt *mdtout,
              mod_mdt_bin_type bin_type);

/** Get an element from an MDT. */
MDTDLLEXPORT
double mdt_get(const struct mod_mdt *mdt, const int indices[], int n_indices,
               GError **err);

/** Set an element in an MDT. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_set(struct mod_mdt *mdt, const int indices[], int n_indices,
                 double val, GError **err);

/** Add mdt2 into mdt1. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_add(struct mdt *mdt1, const struct mdt *mdt2, GError **err);

/** Write input files to plot the given MDT with ASGL. Return TRUE on
    success. */
MDTDLLEXPORT
gboolean mdt_write_asgl(const struct mod_mdt *mdt,
                        const struct mdt_library *mlib, const char *asglroot,
                        const char *text, int dimensions, int every_x_numbered,
                        int every_y_numbered, double plot_density_cutoff,
                        int plots_per_page, int plot_position,
                        const char *plot_type, int x_decimal, int y_decimal,
                        GError **err);

/** Super-duper multi-level hierarchical recursive multi-dimensional
    smoothing of sparse MDT frequency tables. Return TRUE on success. */
MDTDLLEXPORT
gboolean mdt_super_smooth(const struct mdt *mdtin, struct mdt *mdtout,
                          int dimensions, float prior_weight,
                          gboolean entropy_weighing, GError **err);

/** Sum an MDT section. */
MDTDLLEXPORT
double mdt_section_sum(const struct mod_mdt *mdt, const int indices[],
                       int n_indices, GError **err);

/** Get the entropy of an MDT section. */
MDTDLLEXPORT
double mdt_section_entropy(const struct mod_mdt *mdt, const int indices[],
                           int n_indices, GError **err);

/** Get the mean and standard deviation of an MDT section. */
MDTDLLEXPORT
void mdt_section_meanstdev(const struct mod_mdt *mdt,
                           const struct mdt_library *mlib, const int indices[],
                           int n_indices, double *mean, double *stdev,
                           GError **err);

/** Get the version string (static; do not free). */
MDTDLLEXPORT
const char *mdt_version_get(void);

G_END_DECLS

#endif  /* __MDT_MDT_H */
