"""MDT, a module for protein structure analysis.

   MDT prepares a raw frequency table, given information from MODELLER
   alignments and/or PDB files. It can also process the raw frequency table in
   several ways (e.g., L{normalization<mdt.normalize>},
   L{smoothing<mdt.smooth>}), perform "L{entropy<mdt.entropy_full>}"
   calculations, and write out the data in various formats, including
   L{for plotting by ASGL<mdt.write_asgl>} and use as restraints by MODELLER.

   More precisely, MDT uses a sample of sequences, structures, and/or
   alignments to construct a table M{N(a,b,c,...,d)} for features 
   M{a, b, c, ..., d}. The sample for generating the frequencies M{N} is
   obtained depending on the type of features M{a, b, c, ..., d}.
   The sample can contain individual proteins, pairs of proteins, pairs of
   residues in proteins, pairs of aligned residues, pairs of aligned pairs of
   residues, chemical bonds, angles, dihedral angles, and pairs of triplets of
   atoms. Some features work with triple alignments, too. All the needed
   features M{a, b, c, ..., d} are calculated automatically from the sequences,
   alignments, and/or PDB files. The feature bins are defined in a bin file
   which can be changed by the user.

   MDT works by accumulating the table M{N} by processing each sequence or
   alignment in turn. See L{mdt.add_alignment}.

   See U{some studies with MDT<https://salilab.org/internal/manuals/mdt/manual.pdf>} for copious examples.

   @author: Andrej Sali, Ben Webb
   @copyright: 1989-2007 Andrej Sali

   @sort: mdt_library, mdt
"""

__docformat__ = "epytext en"

import _modeller
import _mdt
from modeller.util.modobject import modobject
from modeller.util import modlist

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

    def __init__(self, env, file, binfile, residue_grouping=1,
                 deltai=1, deltaj=1, deltai_ali=False, deltaj_ali=False,
                 distance_atoms=('CA', 'CA'), special_atoms=False,
                 hbond_cutoff=3.5):
        """Create a new MDT library.
           @param env: the Modeller environment to use
           @param file: the MDT definition file
           @param binfile: file defining bin ranges
        """
        self.env = env.copy()
        _modeller.read_mdt_library(self.basept, file)
        _mdt.mdt_library_deltai_set(self.modpt, deltai)
        _mdt.mdt_library_deltaj_set(self.modpt, deltaj)
        _mdt.mdt_library_deltai_ali_set(self.modpt, deltai_ali)
        _mdt.mdt_library_deltaj_ali_set(self.modpt, deltaj_ali)
        _mdt.mdt_library_hbond_cutoff_set(self.modpt, hbond_cutoff)
        _mdt.mdt_library_special_atoms_set(self.modpt, special_atoms)
        _modeller.readbin_mdt_library(self.basept, self.env.libs.modpt,
                                      binfile, residue_grouping, distance_atoms)

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
    def __get_triplet_classes(self):
        return triplet_classes(self)
    def __get_hbond_classes(self):
        return hbond_classes(self)

    modpt = property(__get_modpt)
    basept = property(__get_basept)
    atom_classes = property(__get_atom_classes, doc="Atom classes")
    bond_classes = property(__get_bond_classes, doc="Bond classes")
    angle_classes = property(__get_angle_classes, doc="Angle classes")
    dihedral_classes = property(__get_dihedral_classes, doc="Dihedral classes")
    triplet_classes = property(__get_triplet_classes,
                               doc="Atom triplet classes")
    hbond_classes = property(__get_hbond_classes,
                               doc="Hydrogen bond atom classes")
    deltai = property(__get_deltai, doc="delta i for some feature types")
    deltaj = property(__get_deltaj, doc="delta j for some feature types")
    deltai_ali = property(__get_deltai_ali,
                          doc="True if L{deltai} refers to alignment " + \
                              "positions, or False if residue positions")
    deltaj_ali = property(__get_deltaj_ali,
                          doc="True if L{deltaj} refers to alignment " + \
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



class triplet_classes(bond_classes):
    """Classifications of triplets of atoms into classes"""

    def __init__(self, mlib):
        bond_classes.__init__(self, mlib, 3)

    def read(self, filename):
        """Read atom triplet information from a file"""
        return _modeller.readtriplet_mdt_library(self._mlib.basept, filename)


class hbond_classes(bond_classes):
    """Classifications of atoms into hydrogen bond classes"""

    def __init__(self, mlib):
        bond_classes.__init__(self, mlib, 1)

    def read(self, filename):
        """Read hydrogen bond atom class information from a file"""
        return _mdt.mdt_hbond_read(filename, self._mlib.modpt)


class mdt_section(modobject):
    """A section of a multi-dimensional table"""
    _indices = None
    __mdt = None
    _mlib = None
    _modpt = None

    def __init__(self, mdt, indices):
        self.__mdt = mdt  # Keep a reference to the MDT
        self._modpt = mdt._modpt
        self._mlib = mdt._mlib
        self._indices = indices

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


class mdt(mdt_section):
    """A multi-dimensional table.
       Individual elements from the table can be accessed in standard Python
       fashion, e.g.

       >>> import mdt
       >>> import modeller
       >>> env = modeller.environ()
       >>> mlib = mdt.mdt_library(env, '${LIB}/mdt.lib', '${LIB}/mdt.bin')
       >>> m = mdt.mdt(mlib, features=(1,2,5))
       >>> print m[0,0,0]
    """
    _modpt = None
    _mlib = None

    def __new__(cls, *args, **vars):
        obj = mdt_section.__new__(cls)
        obj._modpt = _modeller.new_mdt_type()
        return obj

    def __init__(self, mlib, file=None, features=None):
        """Create a new MDT.
           @param mlib: the MDT library to use
           @param file: if specified, the filename to read the initial table
                        from
           @param features: if specified, a list of feature types to initialize
                            the table with
        """
        self._indices = ()
        self._mlib = mlib
        if file:
            self.read(file)
        elif features:
            self.make(features)

    def read(self, file):
        """Read an MDT from C{file}."""
        _modeller.read_mdt(self._modpt, self._mlib.basept, file)

    def copy(self):
        """@return: a copy of this MDT.
           @rtype: L{mdt}"""
        mdtout = mdt(self._mlib)
        _modeller.copy_mdt(self._modpt, mdtout._modpt)
        return mdtout

    def make(self, features):
        """Clear the MDT, and set the features"""
        _modeller.make_mdt(self._modpt, self._mlib.basept, features)

    def write(self, file, write_preamble=True):
        """Write an MDT to C{file}. If C{write_preamble} is False, it will
           only write out the contents of the MDT table, without the preamble
           including the feature list, bins, etc. This is useful for example
           for creating a file to be read by another program, such as
           Mathematica."""
        _mdt.mdt_write(self._modpt, self._mlib.modpt, file, write_preamble)

    def reshape(self, features, offset, shape):
        """Reorder the MDT features and optionally decrease their ranges.
           @param features: the new ordering of the MDT features.
           @param offset: the new offset (see L{offset}).
           @param shape: the new shape (see L{shape}).
           @return: the reshaped MDT.
           @rtype: L{mdt}"""
        mdtout = mdt(self._mlib)
        _mdt.mdt_reshape(self._modpt, mdtout._modpt, features, offset, shape)
        return mdtout

    def smooth(self, dimensions, weight):
        """Smooth the MDT with a uniform prior.
           @return: the smoothed MDT.
           @rtype: L{mdt}"""
        mdtout = mdt(self._mlib)
        _mdt.mdt_smooth(self._modpt, mdtout._modpt, dimensions, weight)
        return mdtout

    def normalize(self, dimensions, dx_dy, to_zero, to_pdf):
        """Normalize or scale the MDT. It does not really matter what the
           contents of the input MDT are; sensible contents includes the raw
           or normalized frequencies.
           @param dimensions: specifies whether a 1D or a 2D table is
                  normalized. More precisely, the input distributions are
                  M{p(x/a, b, c, ...)} if C{dimensions}=1, or
                  M{p(x, y/a, b, c, ...)} if C{dimensions}=2, where y and x are
                  the second to last and last features in the list of features.
           @param dx_dy: widths of the bins (either one or two numbers,
                  depending on C{dimensions}). If the value of either dx or dy
                  is -999, the corresponding bin width is extracted from the MDT
                  data structure (not available for all features).
           @param to_zero: if the histogram is empty, setting this True will set
                  the bin values to zero, and False will yield a uniform
                  distribution. It has no effect when the histogram is not
                  empty.
           @return: the normalized MDT.
           @rtype: L{mdt}"""
        mdtout = mdt(self._mlib)
        _mdt.mdt_normalize(self._modpt, mdtout._modpt, self._mlib.modpt,
                           dimensions, dx_dy, to_zero, to_pdf)
        return mdtout

    def integrate(self, features):
        """Integrate the MDT, and reorder the features. This is useful for
           squeezing large MDT arrays into smaller ones, and also for
           eliminating unwanted features (such as X-ray resolution) in
           preparation for L{mdt.write}.
           @param features: the new features (all must be present in the
                            original MDT).
           @return: the integrated MDT.
           @rtype: L{mdt}"""
        mdtout = mdt(self._mlib)
        _mdt.mdt_integrate(self._modpt, mdtout._modpt, features)
        return mdtout

    def exp_transform(self, offset, expoffset, multiplier, power):
        """Apply an exponential transform to the MDT.
           @return: the transformed MDT.
           @rtype: L{mdt}"""
        mdtout = self.copy()
        _mdt.mdt_exp_transform(mdtout._modpt, offset, expoffset, multiplier,
                               power)
        return mdtout

    def log_transform(self, offset, multiplier, undefined=0.):
        """Apply a log transform to the MDT.
           @return: the transformed MDT.
           @rtype: L{mdt}"""
        mdtout = self.copy()
        _mdt.mdt_log_transform(mdtout._modpt, offset, multiplier, undefined)
        return mdtout

    def linear_transform(self, offset, multiplier):
        """Apply a linear transform to the MDT.
           Each element in the new MDT, M{b}, is obtained from the original
           MDT element M{a}, using the following relation:
           M{b = offset + a * multiplier}.
           @return: the transformed MDT.
           @rtype: L{mdt}"""
        mdtout = self.copy()
        _mdt.mdt_linear_transform(mdtout._modpt, offset, multiplier)
        return mdtout

    def inverse_transform(self, offset, multiplier, undefined=0.):
        """Apply an inverse transform to the MDT.
           @return: the transformed MDT.
           @rtype: L{mdt}"""
        mdtout = self.copy()
        _mdt.mdt_inverse_transform(mdtout._modpt, offset, multiplier, undefined)
        return mdtout

    def offset_min(self, dimensions):
        """Offset the MDT by the minimum value.
           @return: the transformed MDT.
           @rtype: L{mdt}"""
        mdtout = self.copy()
        _mdt.mdt_offset_min(mdtout._modpt, dimensions)
        return mdtout

    def close(self, dimensions):
        """Attempt to 'close' the MDT, so that it is useful for creating splines
           of periodic features.
           @return: the closed MDT.
           @rtype: L{mdt}"""
        mdtout = self.copy()
        _mdt.mdt_close(mdtout._modpt, dimensions)
        return mdtout

    def entropy_full(self):
        """Print full entropy information."""
        return _mdt.mdt_entropy_full(self._modpt, self._mlib.modpt)

    def entropy_hx(self):
        """@return: the entropy of the last dependent variable.

           The MDT is integrated to get a 1D histogram, then normalized by
           the sum of the bin values.
           """
        return _mdt.mdt_entropy_hx(self._modpt)

    def super_smooth(self, prior_weight, entropy_weighing):
        """Multi-level smoothing. This super-smoothes the raw frequencies in
           the MDT using the hierarchical smoothing procedure for 1D histograms
           described in Sali and Blundell, JMB 1993. It was also employed in
           Sali and Overington, Prot Sci. 1994.

           Briefly, the idea is to recursively construct the best possible
           prior distribution for smoothing 1D data M{p(x/a, b, c, ...)}.
           The best prior is a weighted sum (weights optionally based on
           entropy) of the best possible estimate of M{p(x/a, b, ...)}
           integrated over c for each c. Each one of these can itself be
           obtained from a prior and the data, and so on recursively.
           @return: the smoothed MDT.
           @rtype: L{mdt}"""
        mdtout = mdt(self._mlib)
        _mdt.mdt_super_smooth(self._modpt, mdtout._modpt, prior_weight,
                              entropy_weighing)
        return mdtout

    def write_asgl(self, asglroot, text, dimensions, plot_position,
                   plots_per_page, plot_density_cutoff=-1., plot_type='HIST2D',
                   every_x_numbered=1, every_y_numbered=1, x_decimal=1,
                   y_decimal=1):
        """Make input files for ASGL.
           @param asglroot: filename prefix for ASGL TOP script and data files.
           @param text: ASGL command lines that are written for each plot.
           @param dimensions: whether to make 1D or 2D plots.
           @param plot_position: position of the plot on the page, in
                  ASGL convention.
           @param plots_per_page: number of plots per page.
           @param plot_density_cutoff: the minimal sum of the bin values that
                  each plot has to have before it is actually written out;
                  otherwise it is ignored. This helps to avoid wasting paper
                  on empty plots when the MDT array data are sparse.
           @param plot_type: select 'HIST2D' or 'PLOT2D' when C{dimensions}=2.
           @param every_x_numbered: spacing for labels on the X axis.
           @param every_y_numbered: spacing for labels on the Y axis.
           @param x_decimal: the number of decimal places used to write
                             X feature values.
           @param y_decimal: the number of decimal places used to write
                             Y feature values.
        """
        return _mdt.mdt_write_asgl(self._modpt, self._mlib.modpt, asglroot,
                                   text, dimensions, every_x_numbered,
                                   every_y_numbered, plot_density_cutoff,
                                   plots_per_page, plot_position, plot_type,
                                   x_decimal, y_decimal)


    def add_alignment(self, aln, distngh=6.0, surftyp=1, accessibility_type=8,
                      residue_span_range=(-9999, 2, 2, 99999), pairs=1,
                      triples=1, io=None, edat=None):
        """Add data from a Modeller alignment to this MDT.
           @param aln: Modeller alignment.
           @param distngh: distance below which residues are considered
                  neighbors.
           @param surftyp: 1 for PSA contact area, 2 for surface area.
           @param accessibility_type: PSA accessibility type (1-10).
        """
        if io is None:
            io = self._mlib.env.io
        if edat is None:
            edat = self._mlib.env.edat
        _mdt.mdt_add_alignment(self._modpt, self._mlib.modpt, aln.modpt,
                               distngh, False, surftyp, accessibility_type,
                               residue_span_range, pairs, triples, io.modpt,
                               edat.modpt, self._mlib.env.libs.modpt)

    def __getitem__(self, indx):
        if not isinstance(indx, (list, tuple)):
            indx = (indx,)
        if len(indx) < len(self.features):
            return mdt_section(self, indx)
        else:
            return _mdt.mdt_get(self._modpt, indx)

    def __get_pdf(self):
        return _mdt.mdt_type_pdf_get(self._modpt)
    def __get_n_proteins(self):
        return _mdt.mdt_type_n_proteins_get(self._modpt)
    def __get_n_protein_pairs(self):
        return _mdt.mdt_type_n_protein_pairs_get(self._modpt)
    def __get_sample_size(self):
        return _mdt.mdt_type_sample_size_get(self._modpt)
    def __get_features(self):
        return feature_list(self)
    def __get_offset(self):
        return tuple([f.offset for f in self.features])
    def __get_shape(self):
        return tuple([len(f.bins) for f in self.features])

    pdf = property(__get_pdf, doc="Whether this MDT is a PDF")
    features = property(__get_features, doc="Features in this MDT")
    n_proteins = property(__get_n_proteins, doc="Number of proteins")
    n_protein_pairs = property(__get_n_protein_pairs,
                               doc="Number of protein pairs")
    sample_size = property(__get_sample_size, doc="Number of sample points")
    shape = property(__get_shape, doc="Array shape")
    offset = property(__get_offset, doc="Array offsets")


class feature_list(modlist.fixlist):
    """A list of all features in an MDT"""

    def __init__(self, mdt):
        self.__mdt = mdt
        modlist.fixlist.__init__(self)

    def __len__(self):
        return _mdt.mdt_type_nfeat_get(self.__mdt._modpt)

    def _getfunc(self, indx):
        return feature(self.__mdt, indx)


class feature(object):
    """A single feature in an MDT"""

    def __init__(self, mdt, indx):
        self._mdt = mdt
        self._indx = indx

    def __get_ifeat(self):
        ifeat = _mdt.mdt_type_ifeat_get(self._mdt._modpt)
        return _modeller.f_int1_get(ifeat, self._indx)
    def __get_bins(self):
        return bin_list(self)
    def __get_offset(self):
        offset = _mdt.mdt_type_istart_get(self._mdt._modpt)
        return _modeller.f_int1_get(offset, self._indx) - 1
    def __get_periodic(self):
        return _mdt.mdt_feature_is_periodic(self.ifeat)

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
        nbins = _mdt.mdt_type_nbins_get(self._mdt._modpt)
        return _modeller.f_int1_get(nbins, self.__feature._indx)

    def _getfunc(self, indx):
        return bin(self.__feature, indx)


class bin(object):
    """A single bin in a feature"""

    def __init__(self, feature, indx):
        self.__feature = feature
        self.__indx = indx

    def __get_symb(self):
        return _mdt.mdt_bin_symbol_get(self.modpt)

    def __get_range(self):
        return ( _mdt.mdt_bin_rang1_get(self.modpt),
                 _mdt.mdt_bin_rang2_get(self.modpt) )

    def __get_modpt(self):
        nfeat = self.__feature._indx
        mdt = self.__feature._mdt._modpt
        mlib = self.__feature._mdt._mlib.modpt
        return _mdt.mdt_library_bin_get(mdt, mlib, nfeat, self.__indx)

    modpt = property(__get_modpt)
    symbol = property(__get_symb, doc="Bin symbol")
    range = property(__get_range, doc="Bin range")


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
    """Write out a Modeller 2D spline library file from an MDT.
       @param fh: Python file to write to
       @param mdt: input MDT, which should be a 2D table (e.g. phi/psi features)
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
