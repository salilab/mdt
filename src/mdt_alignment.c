/** \file mdt_alignment.c      Functions to add alignment data to MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <stdlib.h>
#include <assert.h>
#include <glib.h>
#include <math.h>

#include "modeller.h"
#include "mod_error.h"
#include "mdt.h"
#include "util.h"
#include "geometry.h"
#include "mdt_index.h"
#include "mdt_property.h"
#include "mdt_tuples.h"
#include "mdt_feature.h"
#include "luzzatiplot.h"

/** A source of data for an MDT (generally an alignment) */
struct mdt_source {
  int nseqacc;
  struct mdt_properties *prop;
  struct mod_alignment *aln;
  gboolean *acceptd;
  gboolean sympairs, symtriples;
};

/** For update MDT with error only. Identify the bins needed to be updated,
    and calculate the bincounts for each bin. */
static void gaussian_weight_calc(struct mod_mdt_libfeature *libfeat,
                                 const struct mod_mdt_bin *bin, float numofstd,
                                 int **cpos, float std, float **bincounts,
                                 int **pos, float m1, int *numofbins, int i)
{
  const static float sqrt2 = 1.41421356;
  int j, quitloop=0;
  float lb=0, hb=0,hc,lc;

  lc=m1-numofstd*std;
  hc=m1+numofstd*std;
  for (j = 1; j < libfeat->nbins; j++, bin++) {
    if (lc > bin->rang2) {
      continue;
    } else if (lc < bin->rang1) {
      lb=bin->rang1;
    } else if (lc >= bin->rang1 && lc< bin->rang2) {
      lb=lc;
    }
    if (hc > bin->rang2) {
      hb=bin->rang2;
    } else if (hc >= bin->rang1 && hc < bin->rang2) {
      hb=hc;
      quitloop=1;
    }
    if (*(*(cpos+i)+j-1) >0)  {
      *(*(bincounts+i)+*(*(cpos+i)+j-1)) +=
                    0.5*erf((hb-m1)/(sqrt2*std))
                    -0.5*erf((lb-m1)/(sqrt2*std));
    } else {
      *(*(cpos+i)+j-1)=numofbins[i];
      *(*(pos+i)+numofbins[i])=j;
      *(*(bincounts+i)+numofbins[i]) =
                     0.5*erf((hb-m1)/(sqrt2*std))
                    -0.5*erf((lb-m1)/(sqrt2*std));
      numofbins[i]=numofbins[i]+1;
    }
    if (quitloop>0) break;
  }
}

/** Function to update mdt by taking errors into consideration */
static gboolean update_mdt_witherr(gboolean *outrange, int is1, int ip1,
                                   int is2, int ir1, int ir2, int ir1p,
                                   int ir2p, int ia1, int ia1p,
                                   const struct mdt_library *mlib, int ip2,
                                   struct mdt *mdt, int ibnd1, int ibnd1p,
                                   int is3, int ir3, int ir3p,
                                   const struct mod_libraries *libs,
                                   const struct mod_energy_data *edat,
                                   struct mdt_source *source, GError **err,
                                   float errorscale)
{
  int i,j,*witherr,*periodic,*numofbins, **pos, **cpos,ifi;
  float *mean,*std,**bincounts;
  GError *tmperr = NULL;
  int numofstd=5;
  float m1=0;
  double binval;
  int totalbinnum=1;
  int indx,rn;
  double totalcount=1;
  const struct mdt_tuple *tuple1, *tuple2;
  struct mod_structure *s = mod_alignment_structure_get(source->aln, 0);
  int *indf = g_malloc(sizeof(int) * mdt->base.nfeat);

  *outrange = FALSE;
  witherr = g_malloc0(sizeof(int) * mdt->base.nfeat);
  /* The periodic value is the period of the feature if the feature is
     periodic, or 0 if the feature is not periodic. */
  periodic = g_malloc0(sizeof(int) * mdt->base.nfeat);
  mean = g_malloc0(sizeof(float)* mdt->base.nfeat);
  std = g_malloc0(sizeof(float)* mdt->base.nfeat);
  pos= g_malloc(sizeof(int*)* mdt->base.nfeat);
  cpos= g_malloc(sizeof(int*)* mdt->base.nfeat);
  numofbins = g_malloc0(sizeof(int) * mdt->base.nfeat);
  bincounts= g_malloc(sizeof(float*)* mdt->base.nfeat);

  /** Calculate the mean and error for each feature */
  for (i = 0; i < mdt->base.nfeat && !tmperr && *outrange == FALSE; i++) {
    const struct mod_mdt_feature *feat;
    struct mod_mdt_libfeature *libfeat;
    struct mdt_feature *mfeat;

    feat = &mdt->base.features[i];
    ifi = feat->ifeat;
    libfeat = &mlib->base.features[ifi - 1];
    mfeat = &g_array_index(mlib->features, struct mdt_feature, ifi - 1);

    switch (mfeat->type) {
    case MDT_FEATURE_ATOM_PAIR:
      mean[i] = dist0witherr(ia1, ia1p, s, &std[i], errorscale);
      witherr[i]=1;
      periodic[i]=0;
      break;
    case MDT_FEATURE_TUPLE_PAIR:
      switch (libfeat->name[28]) {
      case 'n':
        mean[i] = dist0witherr(ia1, ia1p, s, &std[i], errorscale);
        periodic[i]=0;
        break;
      case '1':
        tuple2 = property_one_tuple(source->aln, is1, source->prop, mlib,
                                    ibnd1p, ia1p, libs);
        mean[i] = angle0witherr(ia1, ia1p, tuple2->iata[0], s, &std[i],
                                errorscale);
        periodic[i]=180;
        break;
      case '2':
        tuple1 = property_one_tuple(source->aln, is1, source->prop, mlib,
                                    ibnd1, ia1, libs);
        mean[i] = angle0witherr(tuple1->iata[0], ia1, ia1p, s, &std[i],
                                errorscale);
        periodic[i]=180;
        break;
      case 'r':
        tuple1 = property_one_tuple(source->aln, is1, source->prop, mlib,
                                    ibnd1, ia1, libs);
        tuple2 = property_one_tuple(source->aln, is1, source->prop, mlib,
                                    ibnd1p, ia1p, libs);
        mean[i] = dihedral0witherr(tuple1->iata[0], ia1, ia1p, tuple2->iata[0],
                                   s, &std[i], errorscale);
        periodic[i]=360;
        break;
      }
      witherr[i]=1;
      break;
    default:
      numofbins[i]=1;
      *(cpos+i)=g_malloc0(sizeof(int));
      *(pos+i)=g_malloc0(sizeof(int));
      *(bincounts+i)=g_malloc0(sizeof(float));
      *(*(pos+i)+0) = my_mdt_index(ifi, source->aln, is1, ip1, is2, ir1,
                                   ir2, ir1p, ir2p, ia1, ia1p, mlib, ip2,
                                   ibnd1, ibnd1p, is3, ir3, ir3p, libs,
                                   edat, source->prop, &tmperr);
      *(*(bincounts+i)+0)=1;
      witherr[i]=0;
      periodic[i]=0;
      break;
    }
  }

  /* Based on the mean and standard deviation calculated, the bins need to
     be updated and the corresponding bin counts are calculated */
  for (i = 0; i < mdt->base.nfeat && !tmperr && *outrange == FALSE; i++) {
    const struct mod_mdt_feature *feat;
    struct mod_mdt_libfeature *libfeat;
    const struct mod_mdt_bin *bin;
    if (witherr[i]==0) continue;

    feat = &mdt->base.features[i];
    ifi = feat->ifeat;
    libfeat = &mlib->base.features[ifi - 1];
    bin = libfeat->bins;
    numofbins[i]=0;
    *(pos+i)=g_malloc0(sizeof(int)* libfeat->nbins);
    *(cpos+i)=g_malloc0(sizeof(int)* libfeat->nbins);
    *(bincounts+i)=g_malloc0(sizeof(float)* libfeat->nbins);
    m1=mean[i];
    gaussian_weight_calc(libfeat, bin, numofstd, cpos, std[i], bincounts,
                         pos, m1, numofbins, i);
    /* update bin counts again if the feature is periodic*/
    if (periodic[i]>0 && (mean[i]-numofstd*std[i])<bin[0].rang1) {
      if (periodic[i]==180) {
        m1=-mean[i];
      } else if (periodic[i]==360) {
        m1= mean[i]+360;
      }
      gaussian_weight_calc(libfeat, libfeat->bins, numofstd, cpos, std[i],
                           bincounts, pos, m1, numofbins, i);
    }
    /* update bin counts the third time if the feature is periodic*/
    if (periodic[i]>0
        && (mean[i]+numofstd*std[i])>bin[libfeat->nbins-2].rang2) {
      if (periodic[i]==180) {
        m1=360-mean[i];
      } else if (periodic[i]==360) {
        m1= mean[i]-360;
      }
      gaussian_weight_calc(libfeat, libfeat->bins, numofstd, cpos, std[i],
                           bincounts, pos, m1, numofbins, i);
    }
  }

  /* Calculate the total number of bins needed to be updated */
  m1=0;
  for (i = 0; i < mdt->base.nfeat && !tmperr && *outrange == FALSE; i++) {
    totalbinnum=totalbinnum*numofbins[i];
  }

  /* Update individual bins using the bin values calculated above */
  for (i=0; i<totalbinnum; i++) {
    rn=i;
    totalcount=1;
    for (j = 0; j < mdt->base.nfeat; j++) {
      indf[j]=*(*(pos+j)+rn%numofbins[j]);
      totalcount=totalcount*(*(*(bincounts+j)+rn%numofbins[j]));
      rn=rn/numofbins[j];
    }
    indx = indmdt(indf, &mdt->base);
    binval = mod_mdt_bin_get(&mdt->base, indx);
    binval += totalcount;
    m1+=totalcount;
    mod_mdt_bin_set(&mdt->base, indx, binval);
  }

  mdt->sample_size += 1.0;
  for (i=0;i < mdt->base.nfeat; i++) {
    g_free(*(bincounts+i));
    g_free(*(pos+i));
    g_free(*(cpos+i));
  }
  g_free(indf);
  g_free(witherr);
  g_free(periodic);
  g_free(mean);
  g_free(std);
  g_free(pos);
  g_free(cpos);
  g_free(numofbins);
  g_free(bincounts);

  if (tmperr) {
    g_propagate_error(err, tmperr);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Get all MDT indices */
static int *mdt_indices(gboolean *outrange, int is1, int ip1, int is2,
                        int ir1, int ir2, int ir1p, int ir2p, int ia1,
                        int ia1p, const struct mdt_library *mlib, int ip2,
                        const struct mdt *mdt, int ibnd1, int ibnd1p,
                        int is3, int ir3, int ir3p,
                        const struct mod_libraries *libs,
                        const struct mod_energy_data *edat,
                        struct mdt_source *source, GError **err)
{
  int i, *indf;
  GError *tmperr = NULL;
  *outrange = FALSE;
  indf = g_malloc(sizeof(int) * mdt->base.nfeat);

  for (i = 0; i < mdt->base.nfeat && !tmperr && *outrange == FALSE; i++) {
    const struct mod_mdt_feature *feat = &mdt->base.features[i];
    int ifi = feat->ifeat;
    indf[i] = my_mdt_index(ifi, source->aln, is1, ip1, is2, ir1, ir2, ir1p,
                           ir2p, ia1, ia1p, mlib, ip2, ibnd1, ibnd1p, is3,
                           ir3, ir3p, libs, edat, source->prop, &tmperr);
    if (!tmperr && (indf[i] < feat->istart || indf[i] > feat->iend)) {
      *outrange = TRUE;
    }
  }

  if (tmperr) {
    g_propagate_error(err, tmperr);
    return NULL;
  } else {
    return indf;
  }
}

/** Callback function to be called for every valid data point scanned for
    an MDT source. */
typedef gboolean (*mdt_scan_cb) (void *data, struct mdt * mdt, int indx,
                                 GError **err);

/** Simple scan function which updates the counts in the MDT. */
static gboolean scan_update(void *data, struct mdt *mdt, int indx,
                            GError **err)
{
  double binval = mod_mdt_bin_get(&mdt->base, indx);
  binval += 1.0;
  mod_mdt_bin_set(&mdt->base, indx, binval);
  mdt->sample_size += 1.0;
  return TRUE;
}

/** Dummy scan function; used to signal update with errors */
static gboolean scan_update_witherr(void *data, struct mdt *mdt, int indx,
                                    GError **err)
{
  return TRUE;
}

/** Scan function which sums the MDT values. */
static gboolean scan_sum(void *data, struct mdt *mdt, int indx, GError **err)
{
  double *sum = (double *)data;
  *sum += mod_mdt_bin_get(&mdt->base, indx);
  return TRUE;
}

/** Update an MDT with feature data. */
static gboolean update_mdt(struct mdt *mdt, const struct mdt_library *mlib,
                           int is1, int ip1, int is2, int ir1, int ir2,
                           int ir1p, int ir2p, int ip2, int ia1, int ia1p,
                           int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                           const struct mod_libraries *libs,
                           const struct mod_energy_data *edat,
                           struct mdt_source *source, mdt_scan_cb scanfunc,
                           void *scandata, GError **err)
{
  static const char *routine = "update_mdt";
  gboolean outrange;
  int imda, *indf;

  /* Special-casing for updating with errors */
  if (scanfunc == scan_update_witherr) {
    float *errdata = (float*)scandata;
    if (*errdata > 0) {
      return update_mdt_witherr(&outrange, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                                ia1, ia1p, mlib, ip2, mdt, ibnd1, ibnd1p, is3,
                                ir3, ir3p, libs, edat, source, err, *errdata);
    } else {
      /* Fall back to regular function */
      scanfunc = scan_update;
    }
  }

  /* obtain the indices for the feature values in this routine call: */
  indf = mdt_indices(&outrange, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1,
                     ia1p, mlib, ip2, mdt, ibnd1, ibnd1p, is3, ir3, ir3p, libs,
                     edat, source, err);
  if (!indf) {
    return FALSE;
  }

  /* Ignore if any of the indices properly out of range: */
  if (outrange) {
    g_free(indf);
    return TRUE;
  }

  /* obtain the element index for the mdt vector: */
  imda = indmdt(indf, &mdt->base);
  g_free(indf);
  if (imda < 0 || imda >= mdt->base.nelems) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                "%s: MDT index is out of range: %d %d", routine, imda,
                mdt->base.nelems);
    return FALSE;
  }

  return scanfunc(scandata, mdt, imda, err);
}

/** Check to make sure that all features are single protein features */
static gboolean check_single_protein_features(const struct mdt *mdt,
                                              const struct mdt_library *mlib,
                                              GError **err)
{
  int i;
  for (i = 0; i < mdt->base.nfeat; i++) {
    struct mod_mdt_feature *feat = &mdt->base.features[i];
    struct mod_mdt_libfeature *libfeat;
    libfeat = &mlib->base.features[feat->ifeat - 1];
    if (libfeat->iknown > 1) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                  "This combination of features requires a single protein "
                  "scan, but feature number %d (%d) requires multiple "
                  "proteins", i + 1, feat->ifeat);
      return FALSE;
    }
  }
  return TRUE;
}

/** Update the number of protein pairs in the MDT. */
static void update_protein_pairs(struct mdt *mdt, int nseqacc,
                                 gboolean sympairs, gboolean symtriples)
{
  switch (mdt->base.nprotcmp) {
  case 1:
    mdt->n_protein_pairs += nseqacc;
    break;

  case 2:
    if (sympairs) {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) / 2;
    } else {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1);
    }
    break;

  case 3:
    if (symtriples) {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) * (nseqacc - 2) / 6;
    } else {
      mdt->n_protein_pairs += nseqacc * (nseqacc - 1) * (nseqacc - 2);
    }
    break;
  }
}


/** Get start sequence for multiple protein features. */
static int isbeg(int is, gboolean symmetric)
{
  if (symmetric) {
    return is + 1;
  } else {
    return 0;
  }
}

/** Update MDT data for a single protein property. */
static gboolean update_single(struct mdt *mdt, const struct mdt_library *mlib,
                              int is1, int ip1, int ip2, int ir1, int ir1p,
                              const struct mod_libraries *libs,
                              const struct mod_energy_data *edat,
                              struct mdt_source *source,
                              mdt_scan_cb scanfunc, void *scandata,
                              GError **err)
{
  int is2, ir2, ir2p, ia1, ia1p, is3, ir3, ir3p;

  is2 = is3 = is1;
  ir2 = ir3 = ir1;
  ir2p = ir3p = ir1p;
  ia1 = ia1p = 0;

  switch (mdt->scantype) {
    /* whole protein properties tabulated: */
  case MOD_MDTS_PROTEIN:
    if (!update_mdt(mdt, mlib, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ip2,
                    ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, source,
                    scanfunc, scandata, err)) {
      return FALSE;
    }
    break;

    /* residue properties tabulated */
  case MOD_MDTS_RESIDUE:
    if (ir1 >= 0) {
      if (!update_mdt(mdt, mlib, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ip2,
                      ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, source,
                      scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* residue rels compared */
  case MOD_MDTS_RESIDUE_PAIR:
    if (ir1 >= 0 && ir1p >= 0) {
      if (!update_mdt(mdt, mlib, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                      ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat, source,
                      scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;
  }
  return TRUE;
}


/** Update MDT data for a multiple protein property. */
static gboolean update_multiple(struct mdt *mdt,
                                const struct mdt_library *mlib, int is1,
                                int ip1, int ip2, int ir1, int ir1p,
                                struct mdt_source *source,
                                const struct mod_libraries *libs,
                                const struct mod_energy_data *edat,
                                const gboolean acceptd[],
                                mdt_scan_cb scanfunc, void *scandata,
                                GError **err)
{
  int is2, ir2, ir2p, ia1, ia1p, is3, ir3, ir3p;
  const struct mod_alignment *aln = source->aln;

  ia1 = ia1p = 0;

  /* generate all indices for the protein B: */
  for (is2 = isbeg(is1, source->sympairs); is2 < aln->nseq; is2++) {
    if (acceptd[is2]) {
      /* residue indices in the first and second position for protein 2 */
      if (mdt->scantype != MOD_MDTS_PROTEIN) {
        ir2 = mod_int2_get(&aln->ialn, ip1, is2) - 1;
      } else {
        ir2 = 0;
      }
      if (mdt->scantype == MOD_MDTS_RESIDUE_PAIR) {
        ir2p = mod_int2_get(&aln->ialn, ip2, is2) - 1;
      } else {
        ir2p = 0;
      }
      is3 = is2;
      ir3 = ir2;
      ir3p = ir2p;

      /* do you have to generate a pair or a triplet of proteins */
      if (mdt->base.nprotcmp == 2) {
        /* PAIR OF PROTEINS:
           ignore self-comparisons:
           (the second condition only applies in the GETCSR mode to
           allow the prediction from the sequence of the unknown alone
           even if the MDT table was read in that requires the known) */
        if (is1 != is2 || aln->nseq == 1) {
          if (!update_mdt(mdt, mlib, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                          ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs, edat,
                          source, scanfunc, scandata, err)) {
            return FALSE;
          }
        }
      } else {
        /* TRIPLET OF PROTEINS:
           generate all indices for the protein C: */
        for (is3 = isbeg(is2, source->symtriples); is3 < aln->nseq; is3++) {
          if (acceptd[is3]) {
            if (mdt->scantype != MOD_MDTS_PROTEIN) {
              ir3 = mod_int2_get(&aln->ialn, ip1, is3) - 1;
            }
            if (mdt->scantype == MOD_MDTS_RESIDUE_PAIR) {
              ir3p = mod_int2_get(&aln->ialn, ip2, is3) - 1;
            }
            if ((is1 != is2 && is1 != is3 && is2 != is3) || aln->nseq == 1) {
              if (!update_mdt(mdt, mlib, is1, ip1, is2, ir1, ir2, ir1p,
                              ir2p, ip2, ia1, ia1p, 1, 1, is3, ir3, ir3p, libs,
                              edat, source, scanfunc, scandata, err)) {
                return FALSE;
              }
            }
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all proteins or protein pairs in the alignment. */
static gboolean genpair(struct mdt *mdt, const struct mdt_library *mlib,
                        int ip1, int ip2, const struct mod_libraries *libs,
                        const struct mod_energy_data *edat,
                        const gboolean acceptd[], struct mdt_source *source,
                        mdt_scan_cb scanfunc, void *scandata, GError **err)
{
  int is1, ir1, ir1p;
  const struct mod_alignment *aln = source->aln;

  /* generate all indices for protein A: */
  for (is1 = 0; is1 < aln->nseq; is1++) {
    if (acceptd[is1]) {

      /* residue index for a residue of protein A in the 1st position: */
      if (mdt->scantype != MOD_MDTS_PROTEIN) {
        ir1 = mod_int2_get(&aln->ialn, ip1, is1) - 1;
      } else {
        ir1 = 0;
      }
      /* residue index for a residue of protein A in the 2nd position:
         (not used if residue relationships are not compared) */
      if (mdt->scantype == MOD_MDTS_RESIDUE_PAIR) {
        ir1p = mod_int2_get(&aln->ialn, ip2, is1) - 1;
      } else {
        ir1p = 0;
      }

      if (mdt->base.nprotcmp == 1) {
        if (!update_single
            (mdt, mlib, is1, ip1, ip2, ir1, ir1p, libs, edat, source,
             scanfunc, scandata, err)) {
          return FALSE;
        }
      } else {
        if (!update_multiple(mdt, mlib, is1, ip1, ip2, ir1, ir1p, source,
                             libs, edat, acceptd, scanfunc, scandata, err)) {
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}


/** Return TRUE iff the two residues satisfy the residue_span_range check */
static gboolean check_sequence_separation(int ir1, int ir1p,
                                          const int rsrang[4])
{
  int nr = ir1p - ir1;
  return ((nr >= rsrang[0] && nr <= rsrang[1])
          || (nr >= rsrang[2] && nr <= rsrang[3]));
}

/** Return TRUE iff the two chains satisfy the chain_span_range check */
static gboolean check_chain_separation(int chain, int chainp,
                                       const int chain_span_range[4])
{
  return check_sequence_separation(chain, chainp, chain_span_range);
}

/** Return TRUE iff the given atom pair is excluded */
static gboolean check_atom_pair_excluded(int ia1, int ia1p,
                                         GHashTable *exclusions)
{
  return g_hash_table_lookup(exclusions, MAKE_HASH_KEY(ia1, ia1p)) != NULL;
}

/** Scan all residue pairs in the alignment sequence(s). */
static gboolean gen_residue_pairs(struct mdt *mdt,
                                  const struct mdt_library *mlib,
                                  const int rsrang[4],
                                  const struct mod_libraries *libs,
                                  const struct mod_energy_data *edat,
                                  const gboolean acceptd[],
                                  struct mdt_source *source,
                                  mdt_scan_cb scanfunc, void *scandata,
                                  GError **err)
{
  int ip1, ip2;
  const struct mod_alignment *aln = source->aln;

  for (ip1 = 0; ip1 < aln->naln; ip1++) {

    /* only if any of the residue relationships is asymmetric, go NxN
       (mainchain H-bonds are an example) */
    if (mdt->symmetric) {
      for (ip2 = ip1 + rsrang[2]; ip2 <= MIN(aln->naln - 1, ip1 + rsrang[3]);
           ip2++) {
        if (!genpair(mdt, mlib, ip1, ip2, libs, edat, acceptd, source,
                     scanfunc, scandata, err)) {
          return FALSE;
        }
      }
    } else {
      for (ip2 = MAX(0, ip1 + rsrang[0]);
           ip2 <= MIN(aln->naln - 1, ip1 + rsrang[3]); ip2++) {
        if (check_sequence_separation(ip1, ip2, rsrang)) {
          if (!genpair(mdt, mlib, ip1, ip2, libs, edat, acceptd, source,
                       scanfunc, scandata, err)) {
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all atoms in the first alignment sequence. */
static gboolean gen_atoms(struct mdt *mdt, const struct mdt_library *mlib,
                          int is1, const struct mod_libraries *libs,
                          const struct mod_energy_data *edat,
                          struct mdt_source *source, mdt_scan_cb scanfunc,
                          void *scandata, GError **err)
{
  int ia1, ir1, *iresatm;
  struct mod_structure *s1;

  s1 = mod_alignment_structure_get(source->aln, is1);

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    if (!update_mdt(mdt, mlib, is1, 1, 1, ir1, 1, 1, 1, 1, ia1, 1, 1, 1,
                    1, 1, 1, libs, edat, source, scanfunc, scandata, err)) {
      return FALSE;
    }
  }
  return TRUE;
}

/** Add MDT data for a single atom pair. Return TRUE on success. */
static gboolean gen_atompair(struct mdt *mdt, const struct mdt_library *mlib,
                             const int rsrang[4], const int chain_span_range[4],
                             int is1, int ir1, int ia1,
                             int ia1p, int chain, const int iresatm[],
                             const struct mod_libraries *libs,
                             const struct mod_energy_data *edat,
                             struct mdt_source *source, mdt_scan_cb scanfunc,
                             void *scandata, const struct mod_sequence *seq,
                             GHashTable *exclusions, GError **err)
{
  int ir1p = iresatm[ia1p] - 1;
  int chainp = mod_sequence_chain_for_res(seq, ir1p);
  if (check_sequence_separation(ir1, ir1p, rsrang)
      && check_chain_separation(chain, chainp, chain_span_range)
      && (!exclusions || !check_atom_pair_excluded(ia1, ia1p, exclusions))) {
    if (!update_mdt(mdt, mlib, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1,
                    ia1p, 1, 1, 1, 1, 1, libs, edat, source, scanfunc,
                    scandata, err)) {
      return FALSE;
    }
  }
  return TRUE;
}

/** Scan all atom pairs in the first alignment sequence. */
static gboolean gen_atom_pairs(struct mdt *mdt, const struct mdt_library *mlib,
                               const int rsrang[4],
                               const int chain_span_range[4],
                               gboolean exclude_bonds,
                               gboolean exclude_angles,
                               gboolean exclude_dihedrals,
                               int is1, const struct mod_libraries *libs,
                               const struct mod_energy_data *edat,
                               struct mdt_source *source,
                               mdt_scan_cb scanfunc, void *scandata,
                               GError **err)
{
  int ia1, ia1p, ir1, *iresatm;
  struct mod_structure *s1;
  struct mod_sequence *seq;
  GHashTable *exclusions;

  s1 = mod_alignment_structure_get(source->aln, is1);
  seq = mod_alignment_sequence_get(source->aln, is1);

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  iresatm = mod_int1_pt(&s1->cd.iresatm);

  if (!property_exclusions(source->aln, is1, source->prop, mlib,
                           exclude_bonds, exclude_angles, exclude_dihedrals,
                           libs, &exclusions, err)) {
    return FALSE;
  }

  if (mdt->symmetric) {
    for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
      int chain;
      ir1 = iresatm[ia1] - 1;
      chain = mod_sequence_chain_for_res(seq, ir1);
      for (ia1p = ia1 + 1; ia1p < s1->cd.natm; ia1p++) {
        if (!gen_atompair(mdt, mlib, rsrang, chain_span_range, is1, ir1, ia1,
                          ia1p, chain, iresatm, libs, edat, source, scanfunc,
                          scandata, seq, exclusions, err)) {
          return FALSE;
        }
      }
    }
  } else {
    for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
      int chain;
      ir1 = iresatm[ia1] - 1;
      chain = mod_sequence_chain_for_res(seq, ir1);
      for (ia1p = 0; ia1p < s1->cd.natm; ia1p++) {
        if (ia1 != ia1p) {
          if (!gen_atompair(mdt, mlib, rsrang, chain_span_range, is1, ir1, ia1,
                            ia1p, chain, iresatm, libs, edat, source, scanfunc,
                            scandata, seq, exclusions, err)) {
            return FALSE;
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all bonds, angles or dihedrals in the first alignment sequence. */
static gboolean gen_bonds(struct mdt *mdt, const struct mdt_library *mlib,
                          int is1, int npnt, const struct mod_libraries *libs,
                          const struct mod_energy_data *edat,
                          struct mdt_source *source, mdt_scan_cb scanfunc,
                          void *scandata, GError **err)
{
  const struct mdt_bond_list *bonds;
  int ibnd1, is2;

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  is2 = is1;
  bonds = property_bonds(source->aln, is1, source->prop, mlib, npnt, libs);
  for (ibnd1 = 0; ibnd1 < bonds->nbonds; ibnd1++) {
    if (!update_mdt(mdt, mlib, is1, 1, is2, 1, 1, 1, 1, 1, 1, 1, ibnd1,
                    1, 1, 1, 1, libs, edat, source, scanfunc, scandata, err)) {
      return FALSE;
    }
  }
  return TRUE;
}


/** Scan all atom tuples in the first alignment sequence. */
static gboolean gen_atom_tuples(struct mdt *mdt,
                                const struct mdt_library *mlib, int is1,
                                const struct mod_libraries *libs,
                                const struct mod_energy_data *edat,
                                struct mdt_source *source,
                                mdt_scan_cb scanfunc, void *scandata,
                                GError **err)
{
  int ia1, ir1, ibnd1, ibnd1p, ia1p, ir1p, *iresatm;
  struct mod_structure *s1;
  const struct mdt_tuple_list *tup;

  if (!check_single_protein_features(mdt, mlib, err)) {
    return FALSE;
  }
  s1 = mod_alignment_structure_get(source->aln, is1);
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  tup = property_tuples(source->aln, is1, source->prop, mlib, libs);
  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    ir1 = iresatm[ia1] - 1;
    for (ibnd1 = 0; ibnd1 < tup[ia1].ntuples; ibnd1++) {
      /* Just in case you use a single atom feature at position 2 in
         protein A: */
      ia1p = ia1;
      ibnd1p = ibnd1;
      ir1p = ir1;
      if (!update_mdt(mdt, mlib, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1, ia1p,
                      ibnd1, ibnd1p, 1, 1, 1, libs, edat, source, scanfunc,
                      scandata, err)) {
        return FALSE;
      }
    }
  }
  return TRUE;
}


/** Scan all atom tuple pairs in the first alignment sequence. */
static gboolean gen_atom_tuple_pairs(struct mdt *mdt,
                                     const struct mdt_library *mlib,
                                     const int rsrang[4],
                                     const int chain_span_range[4],
                                     gboolean exclude_bonds,
                                     gboolean exclude_angles,
                                     gboolean exclude_dihedrals, int is1,
                                     const struct mod_libraries *libs,
                                     const struct mod_energy_data *edat,
                                     struct mdt_source *source,
                                     mdt_scan_cb scanfunc, void *scandata,
                                     GError **err)
{
  int ia1, ir1, ibnd1, ibnd1p, ia1p, ir1p, *iresatm;
  struct mod_structure *s1;
  struct mod_sequence *seq;
  const struct mdt_tuple_list *tup;
  GHashTable *exclusions;

  s1 = mod_alignment_structure_get(source->aln, is1);
  seq = mod_alignment_sequence_get(source->aln, is1);
  iresatm = mod_int1_pt(&s1->cd.iresatm);
  tup = property_tuples(source->aln, is1, source->prop, mlib, libs);

  if (!property_exclusions(source->aln, is1, source->prop, mlib,
                           exclude_bonds, exclude_angles, exclude_dihedrals,
                           libs, &exclusions, err)) {
    return FALSE;
  }

  for (ia1 = 0; ia1 < s1->cd.natm; ia1++) {
    int chain;
    ir1 = iresatm[ia1] - 1;
    chain = mod_sequence_chain_for_res(seq, ir1);
    for (ibnd1 = 0; ibnd1 < tup[ia1].ntuples; ibnd1++) {
      for (ia1p = 0; ia1p < s1->cd.natm; ia1p++) {
        int chainp;
        ir1p = iresatm[ia1p] - 1;
        chainp = mod_sequence_chain_for_res(seq, ir1p);

        /* the same conditions on sequence separation as for residue pairs */
        if (ia1 != ia1p && check_sequence_separation(ir1, ir1p, rsrang)
            && check_chain_separation(chain, chainp, chain_span_range)
            && (!exclusions
                || !check_atom_pair_excluded(ia1, ia1p, exclusions))) {
          for (ibnd1p = 0; ibnd1p < tup[ia1p].ntuples; ibnd1p++) {
            if (!update_mdt(mdt, mlib, is1, 1, 1, ir1, 1, ir1p, 1, 1, ia1,
                            ia1p, ibnd1, ibnd1p, 1, 1, 1, libs, edat, source,
                            scanfunc, scandata, err)) {
              return FALSE;
            }
          }
        }
      }
    }
  }
  return TRUE;
}


/** Scan all alignment positions or all alignment position pairs in the
    current alignment. If whole protein features only occur in the current
    MDT, then no positions are scanned. */
static gboolean mdt_source_scan(struct mdt *mdt,
                                const struct mdt_library *mlib,
                                struct mdt_source *source, const int rsrang[4],
                                const int chain_span_range[4],
                                gboolean exclude_bonds,
                                gboolean exclude_angles,
                                gboolean exclude_dihedrals,
                                const struct mod_libraries *libs,
                                const struct mod_energy_data *edat,
                                const gboolean acceptd[], int nseqacc,
                                mdt_scan_cb scanfunc, void *scandata,
                                GError **err)
{
  int ip1;

  if (source->aln->nseq == 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "alignment contains no sequences!");
    return FALSE;
  }

  update_protein_pairs(mdt, nseqacc, source->sympairs, source->symtriples);

  switch (mdt->scantype) {
    /* Whole proteins */
  case MOD_MDTS_PROTEIN:
    if (!genpair(mdt, mlib, 1, 1, libs, edat, acceptd, source, scanfunc,
                 scandata, err)) {
      return FALSE;
    }
    break;

    /* Single residues or selected (one per residue) atoms */
  case MOD_MDTS_RESIDUE:
    for (ip1 = 0; ip1 < source->aln->naln; ip1++) {
      if (!genpair(mdt, mlib, ip1, ip1, libs, edat, acceptd, source, scanfunc,
                   scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* intra-molecular residue or selected (one per residue) atom pairs */
  case MOD_MDTS_RESIDUE_PAIR:
    if (!gen_residue_pairs(mdt, mlib, rsrang, libs, edat, acceptd, source,
                           scanfunc, scandata, err)) {
      return FALSE;
    }
    break;

    /* Single protein, all atoms; using only the first protein in
       an alignment! */
  case MOD_MDTS_ATOM:
    if (acceptd[0]) {
      if (!gen_atoms(mdt, mlib, 0, libs, edat, source, scanfunc, scandata,
                     err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom pairs; using only the first protein in
       an alignment! */
  case MOD_MDTS_ATOM_PAIR:
    if (acceptd[0]) {
      if (!gen_atom_pairs(mdt, mlib, rsrang, chain_span_range,
                          exclude_bonds, exclude_angles, exclude_dihedrals,
                          0, libs, edat, source, scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom tuples; using only the first protein in
       an alignment! */
  case MOD_MDTS_TUPLE:
    if (acceptd[0]) {
      if (!gen_atom_tuples(mdt, mlib, 0, libs, edat, source, scanfunc,
                           scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Single protein, all atom tuple pairs; using only the first protein in
       an alignment! */
  case MOD_MDTS_TUPLE_PAIR:
    if (acceptd[0]) {
      if (!gen_atom_tuple_pairs(mdt, mlib, rsrang, chain_span_range,
                                exclude_bonds, exclude_angles,
                                exclude_dihedrals, 0, libs, edat, source,
                                scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all bonds: */
  case MOD_MDTS_BOND:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, 0, MDT_BOND_TYPE_BOND, libs, edat, source,
                     scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all angles: */
  case MOD_MDTS_ANGLE:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, 0, MDT_BOND_TYPE_ANGLE, libs, edat, source,
                     scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;

    /* Scan over all dihedrals: */
  case MOD_MDTS_DIHEDRAL:
    if (acceptd[0]) {
      if (!gen_bonds(mdt, mlib, 0, MDT_BOND_TYPE_DIHEDRAL, libs, edat, source,
                     scanfunc, scandata, err)) {
        return FALSE;
      }
    }
    break;
  }
  return TRUE;
}


/** Prepare a source alignment to add data to an MDT. Returns a source pointer
    (to be later freed with mdt_alignment_close()), or NULL on error. */
struct mdt_source *mdt_alignment_open(struct mdt *mdt,
                                      const struct mdt_library *mlib,
                                      struct mod_alignment *aln, float distngh,
                                      gboolean sdchngh, int surftyp,
                                      int iacc1typ, gboolean sympairs,
                                      gboolean symtriples,
                                      struct mod_io_data *io,
                                      struct mod_libraries *libs, GError **err)
{
  int nseqacc, ierr;

  mod_mdt_getdata(&mdt->base, &nseqacc, aln, distngh, sdchngh, surftyp,
                  iacc1typ, io, libs, &ierr);
  if (ierr) {
    handle_modeller_error(err);
    return NULL;
  }

  mod_lognote("Pre-calculating");
  mod_mdt_precalc(&mdt->base, &mlib->base, aln, libs, &ierr);
  if (ierr) {
    handle_modeller_error(err);
    return NULL;
  } else {
    struct mdt_source *source = g_malloc(sizeof(struct mdt_source));
    source->nseqacc = nseqacc;
    source->aln = aln;
    source->prop = mdt_properties_new(aln);
    source->acceptd = g_malloc(sizeof(gboolean) * aln->nseq);
    source->sympairs = sympairs;
    source->symtriples = symtriples;

    if (mdt->base.readin[0]) {
      int i;
      for (i = 0; i < aln->nseq; i++) {
        struct mod_structure *s1 = mod_alignment_structure_get(aln, i);
        source->acceptd[i] = s1->accepts;
      }
    } else {
      int i;
      for (i = 0; i < aln->nseq; i++) {
        source->acceptd[i] = TRUE;
      }
    }
    return source;
  }
}


/** Close a source alignment previously opened with mdt_alignment_open(). */
void mdt_alignment_close(struct mdt_source *source)
{
  mdt_properties_free(source->prop, source->aln);
  g_free(source->acceptd);
  g_free(source);
}


/** Scan all data points in the source, and return the sum. */
double mdt_source_sum(struct mdt_source *source, struct mdt *mdt,
                      const struct mdt_library *mlib,
                      const int residue_span_range[4],
                      const int chain_span_range[4],
                      gboolean exclude_bonds, gboolean exclude_angles,
                      gboolean exclude_dihedrals,
                      const struct mod_libraries *libs,
                      const struct mod_energy_data *edat, GError **err)
{
  double sum = 0.;
  mdt_source_scan(mdt, mlib, source, residue_span_range, chain_span_range,
                  exclude_bonds, exclude_angles, exclude_dihedrals, libs,
                  edat, source->acceptd, source->nseqacc, scan_sum, &sum, err);
  return sum;
}


/** Return the bin index (starting at 1) of a single MDT feature, at the
    given position in the source alignment. On failure, 0 is returned. */
int mdt_alignment_index(struct mdt_source *source, int ifeat, int is1, int ip1,
                        int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                        int ia1p, int ip2, int ibnd1, int ibnd1p, int is3,
                        int ir3, int ir3p, const struct mdt_library *mlib,
                        const struct mod_libraries *libs,
                        struct mod_energy_data *edat, GError **err)
{
  int indf;
  indf = my_mdt_index(ifeat, source->aln, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                      ia1, ia1p, mlib, ip2, ibnd1, ibnd1p, is3, ir3,
                      ir3p, libs, edat, source->prop, err);
  return indf;
}


/** Add data from an alignment to an MDT. Return TRUE on success. */
gboolean mdt_add_alignment(struct mdt *mdt, const struct mdt_library *mlib,
                           struct mod_alignment *aln, float distngh,
                           gboolean sdchngh, int surftyp, int iacc1typ,
                           const int residue_span_range[4],
                           const int chain_span_range[4],
                           gboolean exclude_bonds, gboolean exclude_angles,
                           gboolean exclude_dihedrals, gboolean sympairs,
                           gboolean symtriples, struct mod_io_data *io,
                           struct mod_energy_data *edat,
                           struct mod_libraries *libs, GError **err)
{
  struct mdt_source *source;

  mod_lognote("Calculating and checking other data: %d", aln->nseq);

  source = mdt_alignment_open(mdt, mlib, aln, distngh, sdchngh, surftyp,
                              iacc1typ, sympairs, symtriples, io, libs, err);
  if (source) {
    gboolean ret;

    mdt->nalns++;
    mdt->n_proteins += source->nseqacc;

    mod_lognote("Updating the statistics array:");
    ret = mdt_source_scan(mdt, mlib, source, residue_span_range,
                          chain_span_range, exclude_bonds, exclude_angles,
                          exclude_dihedrals, libs, edat, source->acceptd,
                          source->nseqacc, scan_update, NULL, err);

    mdt_alignment_close(source);
    return ret;
  } else {
    return FALSE;
  }
}

/** Add data from an alignment to an MDT, with error. Return TRUE on success. */
gboolean mdt_add_alignment_witherr(struct mdt *mdt,
                                   const struct mdt_library *mlib,
                                   struct mod_alignment *aln, float distngh,
                                   gboolean sdchngh, int surftyp, int iacc1typ,
                                   const int residue_span_range[4],
                                   const int chain_span_range[4],
                                   gboolean exclude_bonds,
                                   gboolean exclude_angles,
                                   gboolean exclude_dihedrals,
                                   gboolean sympairs, gboolean symtriples,
                                   struct mod_io_data *io,
                                   struct mod_energy_data *edat,
                                   struct mod_libraries *libs, GError **err,
                                   float errorscale)
{
  struct mdt_source *source;
  float *biso, bisomax=0, bisomaxerror;
  int i, ind;

  mod_lognote("Calculating and checking other data: %d", aln->nseq);

  source = mdt_alignment_open(mdt, mlib, aln, distngh, sdchngh, surftyp,
                              iacc1typ, sympairs, symtriples, io, libs, err);
  if (source) {
    gboolean ret;

    /* The errorscale, used to scale errors, is calculated by assuming the
       atom with the largest Biso has the error defined by R-factor,
       X-ray resolution and the Luzzati plot. */
    struct mod_sequence *seq = mod_alignment_sequence_get(source->aln, 0);
    struct mod_structure *struc = mod_alignment_structure_get(source->aln, 0);
    biso = mod_float1_pt(&struc->cd.biso);
    for(i=0; i< struc->cd.natm; i++){
      if (biso[i]>bisomax) {
        bisomax=biso[i];
      }
    }
    bisomaxerror=sqrt(bisomax/79);
    ind=(int)floor((seq->rfactr)*1000+0.5);
    errorscale=errorscale*bisomaxerror/(luzzatiplot(ind)*(seq->resol));

    mdt->nalns++;
    mdt->n_proteins += source->nseqacc;

    mod_lognote("Updating the statistics array:");
    ret = mdt_source_scan(mdt, mlib, source, residue_span_range,
                          chain_span_range, exclude_bonds, exclude_angles,
                          exclude_dihedrals, libs, edat, source->acceptd,
                          source->nseqacc, scan_update_witherr, &errorscale,
                          err);
    mdt_alignment_close(source);
    return ret;
  } else {
    return FALSE;
  }
}
