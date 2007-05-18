/** \file mdt.h            Functions to handle MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_MDT_H
#define __MDT_MDT_H

#include <glib.h>
#include "mdt_error.h"
#include "mdt_types.h"
#include "mdt_atom_classes.h"
#include "mod_types.h"

G_BEGIN_DECLS

/** Smooth a histogram or the 2D plot with a uniform prior.
    Return TRUE on success. */
gboolean mdt_smooth(const struct mod_mdt *mdtin, struct mod_mdt *mdtout,
                    int dimensions, float weight, GError **err);

/** Transform an MDT with an exp function */
void mdt_exp_transform(struct mod_mdt *mdt, float offset, float expoffset,
                       float multiplier, float power);

/** Transform an MDT with a linear function */
void mdt_linear_transform(struct mod_mdt *mdt, float offset, float multiplier);

/** Transform an MDT with an inverse function */
void mdt_inverse_transform(struct mod_mdt *mdt, float offset, float multiplier,
                           float undefined);

/** Transform an MDT with a log function */
void mdt_log_transform(struct mod_mdt *mdt, float offset, float multiplier,
                       float undefined);

/** Offset an MDT by the minimum value. Return TRUE on success. */
gboolean mdt_offset_min(struct mod_mdt *mdt, int dimensions, GError **err);

/** Close an MDT so that it is useful for creating periodic splines.
    Return TRUE on success. */
gboolean mdt_close(struct mod_mdt *mdt, int dimensions, GError **err);

/** Get the entropy of the dependent variable. */
float mdt_entropy_hx(const struct mod_mdt *mdt, GError **err);

/** Write chi^2, entropies and dependencies of pdf p(x/y,z,...) to the log.
    Return TRUE on success. */
gboolean mdt_entropy_full(const struct mod_mdt *mdt,
                          const struct mdt_library *mlib, GError **err);

/** Write out an MDT. Return TRUE on success. */
gboolean mdt_write(const struct mod_mdt *mdt, const struct mdt_library *mlib,
                   const char *filename, gboolean write_preamble, GError **err);

/** Normalize an MDT. */
gboolean mdt_normalize(const struct mod_mdt *mdtin, struct mod_mdt *mdtout,
                       const struct mdt_library *mlib, int dimensions,
                       const float dx_dy[], int n_dx_dy, gboolean to_zero,
                       gboolean to_pdf, GError **err);

/** Integrate an MDT. Return TRUE on success. */
gboolean mdt_integrate(const struct mod_mdt *mdtin, struct mod_mdt *mdtout,
                       const int features[], int n_features, GError **err);

/** Reshape an MDT. Return TRUE on success. */
gboolean mdt_reshape(const struct mod_mdt *mdtin, struct mod_mdt *mdtout,
                     const int features[], int n_features, const int offset[],
                     int n_offset, const int shape[], int n_shape,
                     GError **err);

/** Clear the MDT array, and set feature types. Return TRUE on success. */
gboolean mdt_make(struct mod_mdt *mdt, const struct mdt_library *mlib,
                  const int features[], int n_features, GError **err);

/** Make mdtout a copy of mdtin. */
void mdt_copy(const struct mod_mdt *mdtin, struct mod_mdt *mdtout);

/** Get an element from an MDT. */
double mdt_get(const struct mod_mdt *mdt, const int indices[], int n_indices,
               GError **err);

/** Write input files to plot the given MDT with ASGL. Return TRUE on
    success. */
gboolean mdt_write_asgl(const struct mod_mdt *mdt,
                        const struct mdt_library *mlib, const char *asglroot,
                        const char *text, int dimensions, int every_x_numbered,
                        int every_y_numbered, double plot_density_cutoff,
                        int plots_per_page, int plot_position,
                        const char *plot_type, int x_decimal, int y_decimal,
                        GError **err);

/** Super-duper multi-level hierarchical recursive multi-dimensional
    smoothing of sparse MDT frequency tables. */
void mdt_super_smooth(const struct mod_mdt *mdtin, struct mod_mdt *mdtout,
                      float prior_weight, gboolean entropy_weighing);

/** Sum an MDT section. */
double mdt_section_sum(const struct mod_mdt *mdt, const int indices[],
                       int n_indices, GError **err);

/** Get the entropy of an MDT section. */
double mdt_section_entropy(const struct mod_mdt *mdt, const int indices[],
                           int n_indices, GError **err);

/** Get the mean and standard deviation of an MDT section. */
void mdt_section_meanstdev(const struct mod_mdt *mdt,
                           const struct mdt_library *mlib, const int indices[],
                           int n_indices, double *mean, double *stdev,
                           GError **err);

/** Is the given feature type periodic? */
gboolean mdt_feature_is_periodic(int ifeat);

/** Add data from an alignment to an MDT. Return TRUE on success. */
gboolean mdt_add_alignment(struct mod_mdt *mdt, const struct mdt_library *mlib,
                           struct alignment *aln, float distngh,
                           gboolean sdchngh, int surftyp, int iacc1typ,
                           const int residue_span_range[4], int pairs,
                           int triples, struct mod_io_data *io,
                           struct mod_energy_data *edat,
                           struct mod_libraries *libs, GError **err);

G_END_DECLS

#endif  /* __MDT_MDT_H */
