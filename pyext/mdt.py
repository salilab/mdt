import _modeller
import _mdt
from modeller.util.modobject import modobject
from modeller.util import modlist

error = _mdt.error

class mdt_library(modobject):
    """Library data used in the construction and use of MDTs"""
    __modpt = None
    env = None

    def __new__(cls, *args, **vars):
        obj = modobject.__new__(cls)
        obj.__modpt = _modeller.new_mdt_library()
        return obj

    def __init__(self, env, file, binfile, residue_grouping=1,
                 deltai=1, deltaj=1, deltai_ali=False, deltaj_ali=False,
                 distance_atoms=('CA', 'CA'), special_atoms=False,
                 hbond_cutoff=3.5):
        self.env = env.copy()
        _modeller.read_mdt_library(self.modpt, file)
        _modeller.mdt_library_deltai_set(self.modpt, deltai)
        _modeller.mdt_library_deltaj_set(self.modpt, deltaj)
        _modeller.mdt_library_deltai_ali_set(self.modpt, deltai_ali)
        _modeller.mdt_library_deltaj_ali_set(self.modpt, deltaj_ali)
        _modeller.readbin_mdt_library(self.modpt, self.env.libs.modpt,
                                      binfile, residue_grouping, distance_atoms,
                                      special_atoms, hbond_cutoff)

    def __del__(self):
        _modeller.free_mdt_library(self.modpt)

    def __get_modpt(self):
        return self.__modpt
    def __get_deltai(self):
        return _modeller.mdt_library_deltai_get(self.modpt)
    def __get_deltaj(self):
        return _modeller.mdt_library_deltaj_get(self.modpt)
    def __get_deltai_ali(self):
        return _modeller.mdt_library_deltai_ali_get(self.modpt)
    def __get_deltaj_ali(self):
        return _modeller.mdt_library_deltaj_ali_get(self.modpt)
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
                          doc="True if deltai refers to alignment " + \
                              "positions, or False if residue positions")
    deltaj_ali = property(__get_deltaj_ali,
                          doc="True if deltaj refers to alignment " + \
                              "positions, or False if residue positions")


class bond_classes(object):
    """Classifications of bonds/angles/dihedrals into classes"""

    def __init__(self, mlib, n_atom):
        self._mlib = mlib
        self.__n_atom = n_atom

    def read(self, filename):
        """Read bond class information from a file"""
        return _modeller.readclass_mdt_library(self._mlib.modpt, filename,
                                               self.__n_atom)


class triplet_classes(bond_classes):
    """Classifications of triplets of atoms into classes"""

    def __init__(self, mlib):
        bond_classes.__init__(self, mlib, 3)

    def read(self, filename):
        """Read atom triplet information from a file"""
        return _modeller.readtriplet_mdt_library(self._mlib.modpt, filename)


class hbond_classes(bond_classes):
    """Classifications of atoms into hydrogen bond classes"""

    def __init__(self, mlib):
        bond_classes.__init__(self, mlib, 1)

    def read(self, filename):
        """Read hydrogen bond atom class information from a file"""
        return _modeller.readhbond_mdt_library(self._mlib.modpt, filename)


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
        return _modeller.mdt_section_sum(self._modpt, self._indices)

    def entropy(self):
        """Entropy of all points in the table"""
        return _modeller.mdt_section_entropy(self._modpt, self._indices)

    def mean_stdev(self):
        """Mean and standard deviation of the table"""
        return _modeller.mdt_section_meanstdev(self._modpt, self._mlib.modpt,
                                               self._indices)


class mdt(mdt_section):
    """A multi-dimensional table"""
    _modpt = None
    _mlib = None

    def __new__(cls, *args, **vars):
        obj = mdt_section.__new__(cls)
        obj._modpt = _modeller.new_mdt_type()
        return obj

    def __init__(self, mlib, file=None, features=None):
        self._indices = ()
        self._mlib = mlib
        if file:
            self.read(file)
        elif features:
            self.make(features)

    def read(self, file):
        """Read an MDT from a file"""
        _modeller.read_mdt(self._modpt, self._mlib.modpt, file)

    def copy(self):
        """Return a copy of this MDT"""
        mdtout = mdt(self._mlib)
        _modeller.copy_mdt(self._modpt, mdtout._modpt)
        return mdtout

    def make(self, features):
        """Clear the MDT, and set the features"""
        _modeller.make_mdt(self._modpt, self._mlib.modpt, features)

    def write(self, file, write_preamble=True):
        """Write an MDT to a file"""
        _mdt.mdt_write(self._modpt, self._mlib.modpt, file, write_preamble)

    def reshape(self, features, offset, shape):
        """Reorder the MDT features and decrease their ranges."""
        mdtout = mdt(self._mlib)
        _modeller.reshape_mdt(self._modpt, mdtout._modpt, features, offset,
                              shape)
        return mdtout

    def smooth(self, dimensions, weight):
        """Smooth the MDT with a uniform prior."""
        mdtout = mdt(self._mlib)
        _mdt.mdt_smooth(self._modpt, mdtout._modpt, dimensions, weight)
        return mdtout

    def normalize(self, dimensions, dx_dy, to_zero, to_pdf):
        """Normalize the MDT."""
        mdtout = mdt(self._mlib)
        _mdt.mdt_normalize(self._modpt, mdtout._modpt, self._mlib.modpt,
                           dimensions, dx_dy, to_zero, to_pdf)
        return mdtout

    def integrate(self, features):
        """Integrate the MDT, and reorder the features."""
        mdtout = mdt(self._mlib)
        _mdt.mdt_integrate(self._modpt, mdtout._modpt, features)
        return mdtout

    def exp_transform(self, offset, expoffset, multiplier, power):
        """Apply an exponential transform to the MDT."""
        mdtout = self.copy()
        _mdt.mdt_exp_transform(mdtout._modpt, offset, expoffset, multiplier,
                               power)
        return mdtout

    def log_transform(self, offset, multiplier, undefined=0.):
        """Apply a log transform to the MDT."""
        mdtout = self.copy()
        _mdt.mdt_log_transform(mdtout._modpt, offset, multiplier, undefined)
        return mdtout

    def linear_transform(self, offset, multiplier):
        """Apply a linear transform to the MDT."""
        mdtout = self.copy()
        _mdt.mdt_linear_transform(mdtout._modpt, offset, multiplier)
        return mdtout

    def inverse_transform(self, offset, multiplier, undefined=0.):
        """Apply an inverse transform to the MDT."""
        mdtout = self.copy()
        _mdt.mdt_inverse_transform(mdtout._modpt, offset, multiplier, undefined)
        return mdtout

    def offset_min(self, dimensions):
        """Offset the MDT by the minimum value."""
        mdtout = self.copy()
        _mdt.mdt_offset_min(mdtout._modpt, dimensions)
        return mdtout

    def close(self, dimensions):
        """Attempt to 'close' the MDT, so that it is useful for creating splines
           of periodic features."""
        mdtout = self.copy()
        _mdt.mdt_close(mdtout._modpt, dimensions)
        return mdtout

    def entropy_full(self):
        """Print full entropy information."""
        return _mdt.mdt_entropy_full(self._modpt, self._mlib.modpt)

    def entropy_hx(self):
        """Get the entropy of the dependent variable."""
        return _mdt.mdt_entropy_hx(self._modpt)

    def super_smooth(self, prior_weight, entropy_weighing):
        """Multi-level smoothing"""
        mdtout = mdt(self._mlib)
        _modeller.super_smooth_mdt(self._modpt, mdtout._modpt, prior_weight,
                                   entropy_weighing)
        return mdtout

    def write_asgl(self, asglroot, text, dimensions, plot_position,
                   plots_per_page, plot_density_cutoff=-1., plot_type='HIST2D',
                   every_x_numbered=1, every_y_numbered=1, x_decimal=1,
                   y_decimal=1):
        """Make input files for ASGL"""
        return _modeller.write_asgl(self._modpt, self._mlib.modpt, asglroot,
                                    text, dimensions, every_x_numbered,
                                    every_y_numbered, plot_density_cutoff,
                                    plots_per_page, plot_position, plot_type,
                                    x_decimal, y_decimal)


    def add_alignment(self, aln, distngh=6.0, surftyp=1, accessibility_type=8,
                      residue_span_range=(-9999, 2, 2, 99999), pairs=1,
                      triples=1, io=None, edat=None):
        """Add data from an alignment to this MDT."""
        if io is None:
            io = self._mlib.env.io
        if edat is None:
            edat = self._mlib.env.edat
        _modeller.add_alignment_mdt(self._modpt, self._mlib.modpt, aln.modpt,
                                    distngh, surftyp, accessibility_type,
                                    residue_span_range, pairs, triples,
                                    io.modpt, edat.modpt,
                                    self._mlib.env.libs.modpt)

    def __getitem__(self, indx):
        if not isinstance(indx, (list, tuple)):
            indx = (indx,)
        if len(indx) < len(self.features):
            return mdt_section(self, indx)
        else:
            return _modeller.mdt_get(self._modpt, indx)

    def __get_pdf(self):
        return _modeller.mdt_type_pdf_get(self._modpt)
    def __get_n_proteins(self):
        return _modeller.mdt_type_n_proteins_get(self._modpt)
    def __get_n_protein_pairs(self):
        return _modeller.mdt_type_n_protein_pairs_get(self._modpt)
    def __get_sample_size(self):
        return _modeller.mdt_type_sample_size_get(self._modpt)
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
        return _modeller.mdt_type_nfeat_get(self.__mdt._modpt)

    def _getfunc(self, indx):
        return feature(self.__mdt, indx)


class feature(object):
    """A single feature in an MDT"""

    def __init__(self, mdt, indx):
        self._mdt = mdt
        self._indx = indx

    def __get_ifeat(self):
        ifeat = _modeller.mdt_type_ifeat_get(self._mdt._modpt)
        return _modeller.f_int1_get(ifeat, self._indx)
    def __get_bins(self):
        return bin_list(self)
    def __get_offset(self):
        offset = _modeller.mdt_type_istart_get(self._mdt._modpt)
        return _modeller.f_int1_get(offset, self._indx) - 1
    def __get_periodic(self):
        return _modeller.mdt_feature_is_periodic(self.ifeat)

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
        nbins = _modeller.mdt_type_nbins_get(self._mdt._modpt)
        return _modeller.f_int1_get(nbins, self.__feature._indx)

    def _getfunc(self, indx):
        return bin(self.__feature, indx)


class bin(object):
    """A single bin in a feature"""

    def __init__(self, feature, indx):
        self.__feature = feature
        self.__indx = indx

    def __get_symb(self):
        nfeat = self.__feature._indx
        mdt = self.__feature._mdt._modpt
        mlib = self.__feature._mdt._mlib.modpt
        return _modeller.mdt_symb_get(mdt, mlib, nfeat, self.__indx)

    def __get_range(self):
        nfeat = self.__feature._indx
        mdt = self.__feature._mdt._modpt
        mlib = self.__feature._mdt._mlib.modpt
        return _modeller.mdt_range_get(mdt, mlib, nfeat, self.__indx)

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
