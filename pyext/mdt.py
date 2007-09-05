"""
   MDT, a module for protein structure analysis.

   MDT prepares a raw frequency table, given information from MODELLER
   alignments and/or PDB files. It can also process the raw frequency table in
   several ways (e.g., normalization with `mdt.normalize`,
   smoothing with `mdt.smooth`), perform entropy calculations with
   `mdt.entropy_full`, and write out the data in various formats, including
   for plotting by ASGL (`mdt.write_asgl`) and use as restraints by MODELLER.

   More precisely, MDT uses a sample of sequences, structures, and/or
   alignments to construct a table *N(a,b,c,...,d)* for features 
   *a, b, c, ..., d*. The sample for generating the frequencies *N* is
   obtained depending on the type of features *a, b, c, ..., d*.
   The sample can contain individual proteins, pairs of proteins, pairs of
   residues in proteins, pairs of aligned residues, pairs of aligned pairs of
   residues, chemical bonds, angles, dihedral angles, and pairs of tuples of
   atoms. Some features work with triple alignments, too. All the needed
   features *a, b, c, ..., d* are calculated automatically from the sequences,
   alignments, and/or PDB files. The feature bins are defined in a bin file
   which can be changed by the user.

   MDT works by accumulating the table *N* by processing each sequence or
   alignment in turn. See `mdt.add_alignment`.

   See `some studies with MDT <https://salilab.org/internal/manuals/mdt/manual.pdf>`__ for copious examples.

   :author: Andrej Sali, Ben Webb
   :copyright: 1989-2007 Andrej Sali
"""

__docformat__ = "restructuredtext"

import _mdt
from modeller.util.modobject import modobject
from modeller.util import modlist
# Import _modeller after modeller itself, since the latter modifies the search
# path for the former:
import _modeller

#: Generic exception
error = _mdt.error

class mdt_library(modobject):
    """Library data used in the construction and use of MDTs"""
    __modpt = None
    env = None

    def __new__(cls, *args, **vars):
        obj = modobject.__new__(cls)
        obj.__modpt = _mdt.mdt_library_new()
        return obj

    def __init__(self, env, binfile, residue_grouping=1,
                 deltai=1, deltaj=1, deltai_ali=False, deltaj_ali=False,
                 distance_atoms=('CA', 'CA'), special_atoms=False,
                 hbond_cutoff=3.5):
        """
        Create a new MDT library.

        :Parameters:
          - `env`: the Modeller environment to use
          - `binfile`: file defining bin ranges
          - `residue_grouping`: type of residue grouping for residue class
            features, as defined in resgrp.lib (1=mainchain conformation,
            2=hydrophobicity)
          - `deltai`: see `deltai`
          - `deltaj`: see `deltaj`
          - `deltai_ali`: see `deltai_ali`
          - `deltaj_ali`: see `deltaj_ali`
          - `distance_atoms`: the atom types to use for "specified" distance
            features
          - `special_atoms`: whether to treat disulfide and termini atoms
            specially for atom class features
          - `hbond_cutoff`: maximum separation between two H-bonded atoms
        """
        self.env = env.copy()
        _mdt.mdt_library_deltai_set(self.modpt, deltai)
        _mdt.mdt_library_deltaj_set(self.modpt, deltaj)
        _mdt.mdt_library_deltai_ali_set(self.modpt, deltai_ali)
        _mdt.mdt_library_deltaj_ali_set(self.modpt, deltaj_ali)
        _mdt.mdt_library_hbond_cutoff_set(self.modpt, hbond_cutoff)
        _mdt.mdt_library_special_atoms_set(self.modpt, special_atoms)
        _modeller.mod_mdt_library_readbin(self.basept, self.env.libs.modpt,
                                          binfile, residue_grouping,
                                          distance_atoms)

    def __del__(self):
        _mdt.mdt_library_free(self.modpt)

    def __get_modpt(self):
        return self.__modpt
    def __get_basept(self):
        return _mdt.mdt_library_base_get(self.__modpt)
    def __get_deltai(self):
        return _mdt.mdt_library_deltai_get(self.modpt)
    def __get_deltaj(self):
        return _mdt.mdt_library_deltaj_get(self.modpt)
    def __get_deltai_ali(self):
        return _mdt.mdt_library_deltai_ali_get(self.modpt)
    def __get_deltaj_ali(self):
        return _mdt.mdt_library_deltaj_ali_get(self.modpt)
    def __get_atom_classes(self):
        return bond_classes(self, 1)
    def __get_bond_classes(self):
        return bond_classes(self, 2)
    def __get_angle_classes(self):
        return bond_classes(self, 3)
    def __get_dihedral_classes(self):
        return bond_classes(self, 4)
    def __get_tuple_classes(self):
        return tuple_classes(self)
    def __get_hbond_classes(self):
        return hbond_classes(self)

    modpt = property(__get_modpt)
    basept = property(__get_basept)
    atom_classes = property(__get_atom_classes, doc="Atom classes")
    bond_classes = property(__get_bond_classes, doc="Bond classes")
    angle_classes = property(__get_angle_classes, doc="Angle classes")
    dihedral_classes = property(__get_dihedral_classes, doc="Dihedral classes")
    tuple_classes = property(__get_tuple_classes,
                             doc="Atom tuple classes")
    hbond_classes = property(__get_hbond_classes,
                               doc="Hydrogen bond atom classes")
    deltai = property(__get_deltai, doc="delta i for some feature types")
    deltaj = property(__get_deltaj, doc="delta j for some feature types")
    deltai_ali = property(__get_deltai_ali,
                          doc="True if `deltai` refers to alignment " + \
                              "positions, or False if residue positions")
    deltaj_ali = property(__get_deltaj_ali,
                          doc="True if `deltaj` refers to alignment " + \
                              "positions, or False if residue positions")


class bond_classes(object):
    """Classifications of bonds/angles/dihedrals into classes"""

    def __init__(self, mlib, n_atom):
        self._mlib = mlib
        self.__n_atom = n_atom

    def read(self, filename):
        """Read bond class information from a file"""
        return _mdt.mdt_atom_classes_read(filename, self._mlib.modpt, 
                                          self.__n_atom)



class tuple_classes(bond_classes):
    """Classifications of tuples of atoms into classes"""

    def __init__(self, mlib):
        bond_classes.__init__(self, mlib, 0)

    def read(self, filename):
        """Read atom tuple information from a file"""
        return _mdt.mdt_tuple_read(filename, self._mlib.modpt)


class hbond_classes(bond_classes):
    """Classifications of atoms into hydrogen bond classes"""

    def __init__(self, mlib):
        bond_classes.__init__(self, mlib, 1)

    def read(self, filename):
        """Read hydrogen bond atom class information from a file"""
        return _mdt.mdt_hbond_read(filename, self._mlib.modpt)


class mdt_section(modobject):
    """A section of a multi-dimensional table"""
    _indices = ()
    __mdt = None
    _mlib = None
    _modpt = None

    def __init__(self, mdt, indices):
        self.__mdt = mdt  # Keep a reference to the MDT
        self._modpt = mdt._modpt
        self._mlib = mdt._mlib
        self._indices = indices

    def _get_removed_rank(self):
        """Return the number of dimensions removed from the full MDT in this
           section (0 for the full MDT, up to nfeatures - 1)"""
        return len(self._indices)

    def sum(self):
        """Sum of all points in the table"""
        return _mdt.mdt_section_sum(self._modpt, self._indices)

    def entropy(self):
        """Entropy of all points in the table"""
        return _mdt.mdt_section_entropy(self._modpt, self._indices)

    def mean_stdev(self):
        """Mean and standard deviation of the table"""
        return _mdt.mdt_section_meanstdev(self._modpt, self._mlib.modpt,
                                          self._indices)

    def __check_indices(self, indices):
        """Make sure the indices for an MDT section are reasonable"""
        for (feat, indx) in zip(self.features, indices):
            istart, iend = feat.offset, len(feat.bins) + feat.offset - 1
            if indx < 0:
                indx += iend + 1
            if not istart <= indx <= iend:
                raise IndexError, "index (%d) not in range %d<=index<=%d" \
                                  % (indx, istart, iend)

    def __getitem__(self, indx):
        if isinstance(indx, list):
            indx = tuple(indx)
        elif not isinstance(indx, tuple):
            indx = (indx,)
        if len(indx) < len(self.features):
            self.__check_indices(indx)
            return mdt_section(self, self._indices + indx)
        else:
            return _mdt.mdt_get(self._modpt, self._indices + indx)

    def __setitem__(self, indx, val):
        if isinstance(indx, list):
            indx = tuple(indx)
        elif not isinstance(indx, tuple):
            indx = (indx,)
        if len(indx) < len(self.features):
            raise ValueError, "Cannot set sections of MDTs"
        else:
            _mdt.mdt_set(self._modpt, self._indices + indx, val)

    def __get_features(self):
        return feature_list(self)
    def __get_offset(self):
        return tuple([f.offset for f in self.features])
    def __get_shape(self):
        return tuple([len(f.bins) for f in self.features])
    features = property(__get_features, doc="Features in this MDT")
    offset = property(__get_offset, doc="Array offsets")
    shape = property(__get_shape, doc="Array shape")


class mdt(mdt_section):
    """A multi-dimensional table.
       Individual elements from the table can be accessed in standard Python
       fashion, e.g.

       >>> import mdt
       >>> import modeller
       >>> env = modeller.environ()
       >>> mlib = mdt.mdt_library(env, '${LIB}/mdt.bin')
       >>> m = mdt.mdt(mlib, features=(1,2,5))
       >>> print m[0,0,0]
    """
    _modpt = None
    _mlib = None

    def __new__(cls, *args, **vars):
        obj = mdt_section.__new__(cls)
        obj._modpt = _modeller.mod_mdt_new()
        return obj

    def __init__(self, mlib, file=None, features=None):
        """
        Create a new MDT.

        :Parameters:
          - `mlib`: the MDT library to use
          - `file`: if specified, the filename to read the initial table from
            (if the name ends with '.hdf5', `mdt.read_hdf5` is used, otherwise
            `mdt.read`)
          - `features`: if specified, a list of feature types to initialize
            the table with
        """
        self._mlib = mlib
        if file:
            if file.endswith(".hdf5"):
                self.read_hdf5(file)
            else:
                self.read(file)
        elif features:
            self.make(features)

    def __del__(self):
        _modeller.mod_mdt_free(self._modpt)

    def read(self, file):
        """Read an MDT from `file`."""
        _modeller.mod_mdt_read(self._modpt, self._mlib.basept, file)

    def read_hdf5(self, file):
        """Read an MDT in HDF5 format from `file`."""
        _mdt.mdt_read_hdf5(self._modpt, self._mlib.modpt, file)

    def copy(self):
        """
        :return: a copy of this MDT.
        :rtype: `mdt`
        """
        mdtout = mdt(self._mlib)
        _mdt.mdt_copy(self._modpt, mdtout._modpt)
        return mdtout

    def make(self, features):
        """Clear the MDT, and set the features"""
        _mdt.mdt_make(self._modpt, self._mlib.modpt, features)

    def write(self, file, write_preamble=True):
        """Write an MDT to `file`. If `write_preamble` is False, it will
           only write out the contents of the MDT table, without the preamble
           including the feature list, bins, etc. This is useful for example
           for creating a file to be read by another program, such as
           Mathematica."""
        _mdt.mdt_write(self._modpt, self._mlib.modpt, file, write_preamble)

    def write_hdf5(self, file):
        """Write an MDT in HDF5 format to `file`."""
        _mdt.mdt_write_hdf5(self._modpt, self._mlib.modpt, file)

    def reshape(self, features, offset, shape):
        """
        Reorder the MDT features and optionally decrease their ranges.

        :Parameters:
          - `features`: the new ordering of the MDT features.
          - `offset`: the new offset (see `offset`).
          - `shape`: the new shape (see `shape`).
        :return: the reshaped MDT.
        :rtype: `mdt`
        """
        mdtout = mdt(self._mlib)
        _mdt.mdt_reshape(self._modpt, mdtout._modpt, features, offset, shape)
        return mdtout

    def smooth(self, dimensions, weight):
        """
        Smooth the MDT with a uniform prior. The MDT is treated either as a
        histogram (if `dimensions` = 1) or a 2D density (`dimensions` = 2),
        and a uniform distribution is added followed by scaling:

        p\ :sub:`i` = |w1| / n + |w2| |vi| / S

        S = |sum|\ :sub:`i`\ :sup:`n` |vi|

        |w1| = 1 / ( 1 + S / (`weight` * n))

        |w2| = 1 - |w1|

        where *v* is the input MDT array, *n* is the number of bins in the
        histogram, and *p* is the output MDT array, smoothed and normalized.
        `weight` is the number of points per bin in the histogram at which
        the relative weights of the input histogram and the uniform prior
        are equal.

        The sum of the bins in the output MDT array is 1, for each histogram.

        Note that the resulting output MDT array is not necessarily a PDF,
        because the bin widths are not taken into account during scaling.
        That is, the sum of all bin values multiplied by the bin widths is not
        1 if the bin widths are not 1.

        :return: the smoothed MDT.
        :rtype: `mdt`

        .. |sum| unicode:: U+03A3
        .. |w1| replace:: w\ :sub:`1`
        .. |w2| replace:: w\ :sub:`2`
        .. |vi| replace:: v\ :sub:`i`
        """
        mdtout = mdt(self._mlib)
        _mdt.mdt_smooth(self._modpt, mdtout._modpt, dimensions, weight)
        return mdtout

    def normalize(self, dimensions, dx_dy, to_zero, to_pdf):
        """
        Normalize or scale the MDT. It does not really matter what the
        contents of the input MDT are; sensible contents include the raw
        or normalized frequencies.

        :Parameters:
          - `dimensions`: specifies whether a 1D or a 2D table is
            normalized. More precisely, the input distributions are
            *p(x|a, b, c, ...)* if `dimensions` = 1, or
            *p(x, y|a, b, c, ...)* if `dimensions` = 2, where y and x are
            the second to last and last features in the list of features.
          - `dx_dy`: widths of the bins (either one or two numbers,
            depending on `dimensions`). If the value of either dx or dy
            is -999, the corresponding bin width is extracted from the MDT
            data structure (not available for all features).
          - `to_zero`: if the histogram is empty, setting this True will set
            the bin values to zero, and False will yield a uniform
            distribution. It has no effect when the histogram is not empty.
          - `to_pdf`: if False, the output is obtained by scaling the input
            such that for 1D histograms |sum| :sub:`i` p(x :sub:`i`) = 1,
            and for 2D histograms |sum| :sub:`i,j` p(x :sub:`i,j`) = 1. Note
            that `dx_dy` is **not** taken into account during this scaling.

            If it is True, the normalization takes into account `dx_dy` so
            that the normalized distribution is actually a PDF. That is,
            |sum| :sub:`i` p(x :sub:`i`) dx = 1 for 1D and 
            |sum| :sub:`i,j` p(x :sub:`i,j`) dx dy = 1 for 2D, where dx and
            dy are the widths of the bins. 
        :return: the normalized MDT.
        :rtype: `mdt`

        .. |sum| unicode:: U+03A3
        """
        mdtout = mdt(self._mlib)
        _mdt.mdt_normalize(self._modpt, mdtout._modpt, self._mlib.modpt,
                           dimensions, dx_dy, to_zero, to_pdf)
        return mdtout

    def integrate(self, features):
        """
        Integrate the MDT, and reorder the features. This is useful for
        squeezing large MDT arrays into smaller ones, and also for
        eliminating unwanted features (such as X-ray resolution) in
        preparation for `mdt.write`.

        :Parameters:
          - `features`: the new features (all must be present in the
            original MDT).
        :return: the integrated MDT.
        :rtype: `mdt`
        """
        mdtout = mdt(self._mlib)
        _mdt.mdt_integrate(self._modpt, mdtout._modpt, features)
        return mdtout

    def exp_transform(self, offset, expoffset, multiplier, power):
        """
        Apply an exponential transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = offset + exp(expoffset + multiplier \* a ^ power)*.

        :rtype: `mdt`
        """
        mdtout = self.copy()
        _mdt.mdt_exp_transform(mdtout._modpt, offset, expoffset, multiplier,
                               power)
        return mdtout

    def log_transform(self, offset, multiplier, undefined=0.):
        """
        Apply a log transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = ln(offset + multiplier \* a)*. Where this would involve the
        logarithm of a negative number, *b* is assigned to be `undefined`.

        :return: the transformed MDT.
        :rtype: `mdt`
        """
        mdtout = self.copy()
        _mdt.mdt_log_transform(mdtout._modpt, offset, multiplier, undefined)
        return mdtout

    def linear_transform(self, offset, multiplier):
        """
        Apply a linear transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = offset + a \* multiplier*.

        :return: the transformed MDT.
        :rtype: `mdt`
        """
        mdtout = self.copy()
        _mdt.mdt_linear_transform(mdtout._modpt, offset, multiplier)
        return mdtout

    def inverse_transform(self, offset, multiplier, undefined=0.):
        """
        Apply an inverse transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = offset + multiplier / a*. Where *a* is zero, *b* is
        assigned to be `undefined`.

        :return: the transformed MDT.
        :rtype: `mdt`
        """
        mdtout = self.copy()
        _mdt.mdt_inverse_transform(mdtout._modpt, offset, multiplier, undefined)
        return mdtout

    def offset_min(self, dimensions):
        """
        Offset the MDT by the minimum value.

        :return: the transformed MDT.
        :rtype: `mdt`
        """
        mdtout = self.copy()
        _mdt.mdt_offset_min(mdtout._modpt, dimensions)
        return mdtout

    def close(self, dimensions):
        """
        Attempt to 'close' the MDT, so that it is useful for creating splines
        of periodic features.

        :return: the closed MDT.
        :rtype: `mdt`
        """
        mdtout = self.copy()
        _mdt.mdt_close(mdtout._modpt, dimensions)
        return mdtout

    def entropy_full(self):
        """Print full entropy information."""
        return _mdt.mdt_entropy_full(self._modpt, self._mlib.modpt)

    def entropy_hx(self):
        """
        The MDT is integrated to get a 1D histogram, then normalized by
        the sum of the bin values. Finally, entropy is calculated as
        |sum|\ :sub:`i` -p\ :sub:`i` ln p\ :sub:`i`

        :return: the entropy of the last dependent variable.
        :rtype: float

        .. |sum| unicode:: U+03A3
        """
        return _mdt.mdt_entropy_hx(self._modpt)

    def super_smooth(self, prior_weight, entropy_weighing):
        """
        Multi-level smoothing. This super-smoothes the raw frequencies in
        the MDT using the hierarchical smoothing procedure for 1D histograms
        described in Sali and Blundell, JMB 1993. It was also employed in
        Sali and Overington, Prot Sci. 1994.

        Briefly, the idea is to recursively construct the best possible
        prior distribution for smoothing 1D data *p(x|a, b, c, ...)*.
        The best prior is a weighted sum (weights optionally based on
        entropy) of the best possible estimate of *p(x|a, b, ...)*
        integrated over *c* for each *c*. Each one of these can itself be
        obtained from a prior and the data, and so on recursively.

        :return: the smoothed MDT.
        :rtype: `mdt`
        """
        mdtout = mdt(self._mlib)
        _mdt.mdt_super_smooth(self._modpt, mdtout._modpt, prior_weight,
                              entropy_weighing)
        return mdtout

    def write_asgl(self, asglroot, text, dimensions, plot_position,
                   plots_per_page, plot_density_cutoff=-1., plot_type='HIST2D',
                   every_x_numbered=1, every_y_numbered=1, x_decimal=1,
                   y_decimal=1):
        """
        Make input files for ASGL.

        :Parameters:
          - `asglroot`: filename prefix for ASGL TOP script and data files.
          - `text`: ASGL command lines that are written for each plot.
          - `dimensions`: whether to make 1D or 2D plots.
          - `plot_position`: position of the plot on the page, in
            ASGL convention.
          - `plots_per_page`: number of plots per page.
          - `plot_density_cutoff`: the minimal sum of the bin values that
            each plot has to have before it is actually written out;
            otherwise it is ignored. This helps to avoid wasting paper
            on empty plots when the MDT array data are sparse.
          - `plot_type`: select 'HIST2D' or 'PLOT2D' when `dimensions` = 2.
          - `every_x_numbered`: spacing for labels on the X axis.
          - `every_y_numbered`: spacing for labels on the Y axis.
          - `x_decimal`: the number of decimal places used to write
            X feature values.
          - `y_decimal`: the number of decimal places used to write
            Y feature values.
        """
        return _mdt.mdt_write_asgl(self._modpt, self._mlib.modpt, asglroot,
                                   text, dimensions, every_x_numbered,
                                   every_y_numbered, plot_density_cutoff,
                                   plots_per_page, plot_position, plot_type,
                                   x_decimal, y_decimal)


    def add_alignment(self, aln, distngh=6.0, surftyp=1, accessibility_type=8,
                      residue_span_range=(-99999, -2, 2, 99999), sympairs=False,
                      symtriples=False, io=None, edat=None):
        """
        Add data from a Modeller alignment to this MDT.

        :Parameters:
          - `aln`: Modeller alignment.
          - `distngh`: distance below which residues are considered neighbors.
          - `surftyp`: 1 for PSA contact area, 2 for surface area.
          - `accessibility_type`: PSA accessibility type (1-10).
          - `residue_span_range`: sequence separation (inclusive) for
            residue-residue, tuple-tuple, and 'any atom' features. For
            the two residue indices r1 and r2 in the tuple-tuple and any
            atom cases, or two alignment position indices in the
            residue-residue case, the following must be true:

            *residue_span_range[0] <= (r2 - r1) <= residue_span_range[1]*

            *residue_span_range[2] <= (r2 - r1) <= residue_span_range[3]*

            For symmetric residue-residue features, only one condition
            must be met:

            *residue_span_range[2] <= abs(r2 - r1) <= residue_span_range[3]*
          - `sympairs`: if True, all features involving pairs of proteins
            are symmetric.
          - `symtriples`: if True, all features involving triples of proteins
            are symmetric.
        """
        if io is None:
            io = self._mlib.env.io
        if edat is None:
            edat = self._mlib.env.edat
        _mdt.mdt_add_alignment(self._modpt, self._mlib.modpt, aln.modpt,
                               distngh, False, surftyp, accessibility_type,
                               residue_span_range, sympairs, symtriples,
                               io.modpt, edat.modpt, self._mlib.env.libs.modpt)


    def open_alignment(self, aln, distngh=6.0, surftyp=1, accessibility_type=8,
                       sympairs=False, symtriples=False, io=None, edat=None):
        """
        Open a Modeller alignment to allow MDT indices to be queried
        (see `source`). Arguments are as for `add_alignment`.

        :rtype: `source`
        """
        return source(self, self._mlib, aln, distngh, surftyp,
                      accessibility_type, sympairs, symtriples, io, edat)

    def __get_pdf(self):
        return _mdt.mod_mdt_pdf_get(self._modpt)
    def __get_n_proteins(self):
        return _mdt.mod_mdt_n_proteins_get(self._modpt)
    def __get_n_protein_pairs(self):
        return _mdt.mod_mdt_n_protein_pairs_get(self._modpt)
    def __get_sample_size(self):
        return _mdt.mod_mdt_sample_size_get(self._modpt)

    pdf = property(__get_pdf, doc="Whether this MDT is a PDF")
    n_proteins = property(__get_n_proteins, doc="Number of proteins")
    n_protein_pairs = property(__get_n_protein_pairs,
                               doc="Number of protein pairs")
    sample_size = property(__get_sample_size, doc="Number of sample points")


class feature_list(modlist.fixlist):
    """A list of all features in an MDT"""

    def __init__(self, mdt):
        self.__mdt = mdt
        self.__removed_rank = mdt._get_removed_rank()
        modlist.fixlist.__init__(self)

    def __len__(self):
        return _mdt.mod_mdt_nfeat_get(self.__mdt._modpt) - self.__removed_rank

    def _getfunc(self, indx):
        return feature(self.__mdt, indx + self.__removed_rank)


class feature(object):
    """A single feature in an MDT"""

    def __init__(self, mdt, indx):
        self._mdt = mdt
        self._indx = indx

    def __get_ifeat(self):
        return _mdt.mod_mdt_feature_ifeat_get(self.modpt)
    def __get_bins(self):
        return bin_list(self)
    def __get_offset(self):
        return _mdt.mod_mdt_feature_istart_get(self.modpt) - 1
    def __get_periodic(self):
        return _mdt.mdt_feature_is_periodic(self.ifeat)
    def __get_modpt(self):
        return _mdt.mod_mdt_feature_get(self._mdt._modpt, self._indx)

    modpt = property(__get_modpt)
    ifeat = property(__get_ifeat, doc="Integer type")
    bins = property(__get_bins, doc="Feature bins")
    offset = property(__get_offset, doc="Offset of first bin")
    periodic = property(__get_periodic, doc="Whether feature is periodic")


class bin_list(modlist.fixlist):
    """A list of all bins in a feature"""

    def __init__(self, feature):
        self.__feature = feature
        self._mdt = feature._mdt
        modlist.fixlist.__init__(self)

    def __len__(self):
        return _mdt.mod_mdt_feature_nbins_get(self.__feature.modpt)

    def _getfunc(self, indx):
        return bin(self.__feature, indx)


class bin(object):
    """A single bin in a feature"""

    def __init__(self, feature, indx):
        self.__feature = feature
        self.__indx = indx

    def __get_symb(self):
        return _mdt.mod_mdt_bin_symbol_get(self.modpt)

    def __get_range(self):
        return ( _mdt.mod_mdt_bin_rang1_get(self.modpt),
                 _mdt.mod_mdt_bin_rang2_get(self.modpt) )

    def __get_modpt(self):
        nfeat = self.__feature._indx
        mdt = self.__feature._mdt._modpt
        mlib = self.__feature._mdt._mlib.modpt
        return _mdt.mdt_library_bin_get(mdt, mlib, nfeat, self.__indx)

    modpt = property(__get_modpt)
    symbol = property(__get_symb, doc="Bin symbol")
    range = property(__get_range, doc="Bin range")


class source(object):
    """A source of data for an MDT (generally a Modeller alignment, opened
       with `mdt.open_alignment()`)."""

    def __init__(self, mdt, mlib, aln, distngh, surftyp, accessibility_type,
                 sympairs, symtriples, io, edat):
        self._mdt = mdt
        self._mlib = mlib
        self._aln = aln
        if io is None:
            io = mlib.env.io
        if edat is None:
            edat = mlib.env.edat
        self._edat = edat
        self._modpt = _mdt.mdt_alignment_open(mdt._modpt, mlib.modpt, aln.modpt,
                                              distngh, False, surftyp,
                                              accessibility_type, sympairs,
                                              symtriples, io.modpt,
                                              mlib.env.libs.modpt)

    def __del__(self):
        if hasattr(self, "_modpt"):
            _mdt.mdt_alignment_close(self._modpt)

    def sum(self, residue_span_range=(-99999, -2, 2, 99999)):
        """Scan all data points in the source, and return the sum."""
        f = _mdt.mdt_source_sum
        return f(self._modpt, self._mdt._modpt, self._mlib.modpt,
                 residue_span_range, self._mlib.env.libs.modpt,
                 self._edat.modpt)

    def index(self, ifeat, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1, ia1p,
              ip2, ibnd1, ibnd1p, is3, ir3, ir3p):
        """
        Return the bin index (starting at 1) of a single MDT feature.
        (Arguments ending in 2 and 3 are used for features involving pairs
        or triples of proteins.)

        :Parameters:
          - `ifeat`: MDT feature type.
          - `is1`: index of the sequence within the alignment.
          - `ip1`: position within the sequence (i.e. including gaps).
          - `ir1`: residue index (i.e. not including alignment gaps).
          - `ir1p`: second residue index for residue-residue features.
          - `ia1`: atom index.
          - `ia1p`: second atom index for atom-atom features.
          - `ibnd1`: bond or tuple index.
          - `ibnd1p`: second bond/tuple index for bond-bond or tuple-tuple
            features.
        """
        f = _mdt.mdt_alignment_index
        return f(self._modpt, ifeat, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1,
                 ia1p, ip2, ibnd1, ibnd1p, is3, ir3, ir3p, self._mlib.modpt,
                 self._mlib.env.libs.modpt, self._edat.modpt)


def _pass_cutoffs(mdt, num, bin, density_cutoff, entropy_cutoff):
    if density_cutoff is not None:
        sum = mdt[num].sum()
        if sum < density_cutoff:
            print "Restraint %s skipped: density %.4f below cutoff %.4f" \
                  % (bin.symbol, sum, density_cutoff)
            return False
    if entropy_cutoff is not None:
        entropy = mdt[num].entropy()
        if entropy > entropy_cutoff:
            print "Restraint %s skipped: entropy %.4f above cutoff %.4f" \
                  % (bin.symbol, entropy, entropy_cutoff)
            return False
    return True


def _write_meanstdevlib(fh, mdt, numat, phystype, feattype, convfunc,
                        density_cutoff=None, entropy_cutoff=None):
    print >> fh, "#   residue    atoms        mean    stdev"
    print >> fh, "_params = ["
    for (num,bin) in enumerate(mdt.features[0].bins):
        if _pass_cutoffs(mdt, num, bin, density_cutoff, entropy_cutoff):
            symbols = bin.symbol.split(':')
            res = symbols[0]
            ats = tuple(symbols[1:])
            if len(ats) != numat:
                raise ValueError, "Bin name %s should be res. plus %d atoms" \
                                  % (bin.symbol, numat)
            mean, stdev = mdt[num].mean_stdev()
            print >> fh, "    ( '%s', %s, %.4f, %.4f )," \
                         % (res, str(ats), convfunc(mean), convfunc(stdev))
    print >> fh, """  ]

def make_restraints(atmsel, restraints, num_selected):
    from modeller import forms, physical, features
    for (res, atoms, mean, stdev) in _params:
        for a in atmsel.find_atoms(res, atoms, num_selected):
            r = forms.gaussian(%s, %s(*a), mean,
                               stdev)
            restraints.add(r)""" % (phystype, feattype)

def _noconv(a):
    return a
def _degrees_to_radians(a):
    pi = 3.14159265358979323846
    return a / 180.0 * pi

def write_bondlib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    _write_meanstdevlib(fh, mdt, 2, "physical.bond", "features.distance",
                        _noconv, density_cutoff, entropy_cutoff)

def write_anglelib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    _write_meanstdevlib(fh, mdt, 3, "physical.angle", "features.angle",
                        _degrees_to_radians, density_cutoff, entropy_cutoff)

def write_improperlib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    _write_meanstdevlib(fh, mdt, 4, "physical.improper", "features.dihedral",
                        _degrees_to_radians, density_cutoff, entropy_cutoff)

def _get_splinerange(feat):
    periodic = feat.periodic

#   histogram bin size in real units:
    startrange = feat.bins[0].range
    endrange = feat.bins[-1].range
    dx = startrange[1] - startrange[0]
#   left and right value of the range of the spline:
    x1 = startrange[0] + 0.5 * dx
    if periodic:
        x2 = endrange[1] + 0.5 * dx
    else:
        x2 = endrange[1] - 0.5 * dx
    dx = _degrees_to_radians(dx)
    x1 = _degrees_to_radians(x1)
    x2 = _degrees_to_radians(x2)
    return periodic, dx, x1, x2

def write_splinelib(fh, mdt, dihtype, density_cutoff=None, entropy_cutoff=None):
    (periodic, dx, x1, x2) = _get_splinerange(mdt.features[1])

    print >> fh, "#   residue   spline values"
    print >> fh, "_params = ["
    for (nx,bin) in enumerate(mdt.features[0].bins):
        if _pass_cutoffs(mdt, nx, bin, density_cutoff, entropy_cutoff):
            splinevals = ["%.4f" % mdt[nx,ny] \
                          for ny in range(len(mdt.features[1].bins))]
            print >> fh, "    ( '%s', (%s) )," \
                         % (bin.symbol, ', '.join(splinevals))
    print >> fh, """  ]

def make_restraints(atmsel, restraints, num_selected):
    from modeller import forms, physical, features
    for (res, values) in _params:
        arr = True
        for a in atmsel.find_%s_dihedrals(res, num_selected):
            r = forms.spline(physical.%s_dihedral,
                             features.dihedral(*a), open=%s, low=%.5f,
                             high=%.5f, delta=%.5f, lowderiv=0,
                             highderiv=0, values=values, use_array=arr)
            arr = restraints.add(r)""" % (dihtype, dihtype, not periodic,
                                          x1, x2, dx)


def write_2dsplinelib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    """
    Write out a Modeller 2D spline library file from an MDT.

    :Parameters:
      - `fh`: Python file to write to
      - `mdt`: input MDT, which should be a 2D table (e.g. phi/psi features)
    """
    (yperiodic, dy, y1, y2) = _get_splinerange(mdt.features[1])
    (zperiodic, dz, z1, z2) = _get_splinerange(mdt.features[2])

    print >> fh, "#   residue   spline values"
    print >> fh, "_params = ["
    for (nx,bin) in enumerate(mdt.features[0].bins):
        if _pass_cutoffs(mdt, nx, bin, density_cutoff, entropy_cutoff):
            splinevals = []
            for ny in range(len(mdt.features[1].bins)):
                for nz in range(len(mdt.features[2].bins)):
                    splinevals.append("%.4f" % mdt[nx,ny,nz])
            print >> fh, "    ( '%s', (%s) )," \
                         % (bin.symbol, ', '.join(splinevals))
    print >> fh, """  ]

def make_restraints(atmsel, restraints, num_selected):
    from modeller import forms, physical, features
    for (res, values) in _params:
        arr = True
        for a in atmsel.find_atoms(res, ('-C', 'N', 'CA', 'C',
                                         'N', 'CA', 'C', '+N'), num_selected):
            r = forms.nd_spline(physical.phi_psi_dihedral, values,
                                dimensions=(%d,%d), use_array=arr)
            r.add_dimension(features.dihedral(*a[:4]), open=%s,
                            low=%.5f, high=%.5f, delta=%.5f,
                            lowderiv=0., highderiv=0.)
            r.add_dimension(features.dihedral(*a[4:]), open=%s,
                            low=%.5f, high=%.5f, delta=%.5f,
                            lowderiv=0., highderiv=0.)
            arr = restraints.add(r)""" % (len(mdt.features[1].bins),
                                          len(mdt.features[2].bins),
                                          not yperiodic, y1, y2, dy,
                                          not zperiodic, z1, z2, dz)
