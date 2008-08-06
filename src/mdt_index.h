/** \file mdt_index.h      Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#ifndef __MDT_INDEX_H
#define __MDT_INDEX_H

#include <glib.h>
#include "mdt_config.h"
#include "mod_types.h"
#include "mdt_types.h"

G_BEGIN_DECLS

/** 3 bond types: bond, angle, dihedral */
#define MDT_BOND_TYPE_BOND     0
#define MDT_BOND_TYPE_ANGLE    1
#define MDT_BOND_TYPE_DIHEDRAL 2
#define N_MDT_BOND_TYPES       3

/** Atom tuple (1 or 2 atoms stored for each leading atom) */
struct mdt_tuple {
  /** Atom indices */
  int iata[2];
  /** Tuple class */
  int tupclass;
};

/** List of atom tuples */
struct mdt_tuple_list {
  /** Number of tuples */
  int ntuples;
  /** Tuple data */
  struct mdt_tuple *tuples;
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

/** Precalculated per-sequence properties for calculating MDT indices */
struct mdt_properties {
  /** Lists of bonds */
  struct mdt_bond_list *bonds[N_MDT_BOND_TYPES];
  /** Lists of atom tuples for each atom */
  struct mdt_tuple_list *tuples;
  /** Bin indices for hydrogen bond atom type */
  int *hb_iatta;
  /** Hydrogen bond satisfaction index */
  float *hbpot;
  /** Radius of gyration, or -1 if not yet calculated */
  float radius_gyration;
  /** Bin indices for atom type */
  int *iatta;
  /** Bin indices for atom accessibility */
  int *iatmacc;
  /** Bin indices for fractional atom accessibility */
  int *ifatmacc;
};

/** Register our MDT feature types */
MDTDLLLOCAL
void mdt_register_features(struct mod_mdt_library *mlib);

/** Calculate a single MDT feature index */
MDTDLLLOCAL
int my_mdt_index(int ifi, const struct mod_alignment *aln, int is1, int ip1,
                 int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                 int ia1p, const struct mdt_library *mlib, int ip2,
                 int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                 const struct mod_libraries *libs,
                 const struct mod_energy_data *edat,
                 struct mdt_properties *prop, GError **err);

/** Convert a raw number to the corresponding feature's MDT bin index */
MDTDLLLOCAL
int iclsbin(float x, const struct mod_mdt_libfeature *feat);

G_END_DECLS

#endif  /* __MDT_INDEX_H */
