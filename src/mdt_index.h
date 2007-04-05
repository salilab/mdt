/** \file mdt_index.h      Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#ifndef __MDT_INDEX_H
#define __MDT_INDEX_H

#include <glib.h>
#include "mod_types.h"
#include "mdt_types.h"

/* Allow building with glib < 2.6 */
#ifndef G_GNUC_INTERNAL
#define G_GNUC_INTERNAL
#endif

G_BEGIN_DECLS

/** 3 bond types: bond, angle, dihedral */
#define MDT_BOND_TYPE_BOND     0
#define MDT_BOND_TYPE_ANGLE    1
#define MDT_BOND_TYPE_DIHEDRAL 2
#define N_MDT_BOND_TYPES       3

/** Atom triplet (2 atoms stored for each leading atom) */
struct mdt_triplet {
  /** Atom indices */
  int iata[2];
  /** Triplet class */
  int trpclass;
};

/** List of atom triplets */
struct mdt_triplet_list {
  /** Number of triplets */
  int ntriplets;
  /** Triplet data */
  struct mdt_triplet *triplets;
};

/** A single bond/angle/dihedral in a template structure */
struct mdt_bond {
  /** Indices of all atoms in the bond */
  int iata[4];
  /** Bond type index */
  int bndgrp;
};

/** A list of bonds */
struct mdt_bond_list {
  /** Number of bonds */
  int nbonds;
  /** The bonds */
  struct mdt_bond *bonds;
};

/** Properties for calculating MDT indices */
struct mdt_properties {
  /** Lists of bonds */
  struct mdt_bond_list *bonds[N_MDT_BOND_TYPES];
  /** Lists of atom triplets for each atom */
  struct mdt_triplet_list *triplets;
  /** Bin indices for hydrogen bond atom type */
  int *hb_iatta;
  /** Hydrogen bond satisfaction index */
  float *hbpot;
  /** Bin indices for atom type */
  int *iatta;
  /** Bin indices for atom accessibility */
  int *iatmacc;
  /** Bin indices for fractional atom accessibility */
  int *ifatmacc;
};

/** Make a new mdt_properties structure */
G_GNUC_INTERNAL
struct mdt_properties *mdt_properties_new(const struct alignment *aln);

/** Free an mdt_properties structure */
G_GNUC_INTERNAL
void mdt_properties_free(struct mdt_properties *prop,
                         const struct alignment *aln);

/** Calculate a single MDT feature index */
G_GNUC_INTERNAL
int my_mdt_index(int ifi, const struct alignment *aln, int is1, int ip1,
                 int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                 int ia1p, const struct mdt_library *mlib, int ip2,
                 int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                 const struct libraries *libs,
                 const struct energy_data *edat,
                 struct mdt_properties *prop, GError **err);

/** Get/calculate the list of all bonds for a structure. */
G_GNUC_INTERNAL
const struct mdt_bond_list *property_bonds(const struct alignment *aln, int is,
                                           struct mdt_properties *prop,
                                           const struct mdt_library *mlib,
                                           int bondtype,
                                           const struct libraries *libs);

/** Get/calculate the list of all triplets for a structure. */
G_GNUC_INTERNAL
const struct mdt_triplet_list *property_triplets(const struct alignment *aln,
                                                 int is,
                                                 struct mdt_properties *prop,
                                                 const struct mdt_library *mlib,
                                                 const struct libraries *libs);

G_END_DECLS

#endif  /* __MDT_INDEX_H */
