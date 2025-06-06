# -*- coding=utf-8 -*-

"""
   MDT, a module for protein structure analysis.

   Copyright 1989-2025 Andrej Sali.

   MDT is free software: you can redistribute it and/or modify
   it under the terms of version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with MDT.  If not, see <http://www.gnu.org/licenses/>.
"""

__docformat__ = "restructuredtext"

__all__ = ['MDTError', 'FileFormatError', 'TableSection', 'Table', 'Library',
           'Feature', 'Bin', 'Source', 'BondClasses', 'TupleClasses',
           'HydrogenBondClasses', 'Float', 'Double', 'Int32', 'UnsignedInt32',
           'Int16', 'UnsignedInt16', 'Int8', 'UnsignedInt8',
           'write_2dsplinelib', 'write_anglelib', 'write_bondlib',
           'write_improperlib', 'write_splinelib', 'uniform_bins',
           'write_statpot']

try:
    from modeller.util.modobject import ModObject
except ImportError:
    from modeller.util.modobject import modobject as ModObject
from modeller.util import modlist
# Import _modeller after modeller itself, since the latter modifies the search
# path for the former:
import _modeller
import _mdt

del _modeller

#: Generic MDT exception
MDTError = _mdt.MDTError

#: File format error
FileFormatError = _mdt.FileFormatError


# Get version info
def __get_version_info(version):
    try:
        return tuple([int(x) for x in version.split('.')])
    except ValueError:
        return version


__version__ = version = _mdt.mdt_version_get()
version_info = __get_version_info(version)


def _prepare_bond_span(bond_span):
    """Helper function for bond_span_range"""
    if bond_span is None:
        return (-1, -1)
    else:
        return bond_span


class _BinType(object):
    def __init__(self, bin_type):
        self._bin_type = bin_type


#: Single-precision floating-point bin storage (approximate range 0 to +/-1e38)
Float = _BinType(_mdt.MOD_MDTB_FLOAT)
#: Double-precision floating-point bin storage (approximate range
#: 0 to +/-1e308)
Double = _BinType(_mdt.MOD_MDTB_DOUBLE)
#: 32-bit signed integer bin storage (range -1e31 to 1e31)
Int32 = _BinType(_mdt.MOD_MDTB_INT32)
#: 32-bit unsigned integer bin storage (range 0 to 1e32)
UnsignedInt32 = _BinType(_mdt.MOD_MDTB_UINT32)
#: 16-bit signed integer bin storage (range -1e15 to 1e15)
Int16 = _BinType(_mdt.MOD_MDTB_INT16)
#: 16-bit unsigned integer bin storage (range 0 to 1e16)
UnsignedInt16 = _BinType(_mdt.MOD_MDTB_UINT16)
#: 8-bit signed integer bin storage (range -127 to 128)
Int8 = _BinType(_mdt.MOD_MDTB_INT8)
#: 8-bit unsigned integer bin storage (range 0 to 255)
UnsignedInt8 = _BinType(_mdt.MOD_MDTB_UINT8)


class Library(ModObject):
    """
    Library data used in the construction and use of MDTs.

    :Parameters:
      - `env`: the Modeller environment to use
      - `distance_atoms`: the atom types to use for the
        :class:`features.ResidueDistance` feature
      - `special_atoms`: whether to treat disulfide and termini atoms
        specially for atom class features (see :class:`features.AtomType`)
      - `hbond_cutoff`: maximum separation between two H-bonded atoms
        (see :class:`features.HydrogenBondDonor`)
    """
    _modpt = None
    _env = None

    def __init__(self, env, distance_atoms=('CA', 'CA'), special_atoms=False,
                 hbond_cutoff=3.5):
        self._env = env.copy()
        self._modpt = _mdt.mdt_library_new(self._env.libs.modpt, self)
        _mdt.mdt_library_hbond_cutoff_set(self._modpt, hbond_cutoff)
        _mdt.mdt_library_special_atoms_set(self._modpt, special_atoms)
        _mdt.mdt_library_distance_atoms_set(self._modpt, *distance_atoms)

    def __del__(self):
        _mdt.mdt_library_free(self._modpt)

    def __get_atom_classes(self):
        return BondClasses(self, 1)

    def __get_bond_classes(self):
        return BondClasses(self, 2)

    def __get_angle_classes(self):
        return BondClasses(self, 3)

    def __get_dihedral_classes(self):
        return BondClasses(self, 4)

    def __get_tuple_classes(self):
        return TupleClasses(self)

    def __get_hbond_classes(self):
        return HydrogenBondClasses(self)

    atom_classes = property(__get_atom_classes,
                            doc="Atom classes; see :class:`BondClasses`")
    bond_classes = property(__get_bond_classes,
                            doc="Bond classes; see :class:`BondClasses`")
    angle_classes = property(__get_angle_classes,
                             doc="Angle classes; see :class:`BondClasses`")
    dihedral_classes = property(
        __get_dihedral_classes,
        doc="Dihedral classes; see :class:`BondClasses`")
    tuple_classes = property(
        __get_tuple_classes,
        doc="Atom tuple classes; see :class:`TupleClasses`"
            " and :ref:`tuple_features`")
    hbond_classes = property(__get_hbond_classes,
                             doc="Hydrogen bond atom classes; "
                                 "see :class:`HydrogenBondClasses`")


class BondClasses(object):
    """Classifications of atoms/bonds/angles/dihedrals into classes.
       These classes are used by
       :ref:`atom <atom_features>` and
       :ref:`chemical bond <chemical_bond_features>` features.
       Usually accessed as :attr:`Library.atom_classes`,
       :attr:`Library.bond_classes`, :attr:`Library.angle_classes`, or
       :attr:`Library.dihedral_classes`. (There is no
       need to create your own BondClasses objects.)"""

    def __init__(self, mlib, n_atom):
        self._mlib = mlib
        self.__n_atom = n_atom

    def read(self, filename):
        """Read class information from `filename`.
           This is a text file with a simple format. Each line either
           denotes the start of a new named class, or names a member of the
           last-named class, as a residue name followed by one or more atom
           names. For example, an atom class file might start with::

             ATMGRP 'AC'
               ATOM 'ALA' 'CA'
               ATOM 'ALA' 'C'
               ATOM '*' 'CB'

           Thus, the first atom class is called 'AC' and any CA or C atom in
           an ALA residue, or the CB atom in any residue, will be placed in
           this class.

           Bond class files are similar but use BNDGRP and BOND lines,
           each of which names two atoms::

             BNDGRP 'ALA:C:+N'
               BOND 'ALA' 'C' '+N'

           Note that CHARMM-style + or - prefixes can be added to atom names
           for all but the first atom on a BOND line, to indicate the atom
           must be found in the next or previous residue.

           Angle class files use ANGGRP and ANGLE lines; each ANGLE line
           names three atoms. Dihedral class files use DIHGRP and DIHEDRAL
           lines; each DIHEDRAL line names four atoms.
        """
        return _mdt.mdt_atom_classes_read(filename, self._mlib._modpt,
                                          self.__n_atom)


class TupleClasses(BondClasses):
    """Classifications of tuples of atoms into classes.
       Usually accessed as :attr:`Library.tuple_classes`.
       These classes are used by :ref:`tuple <tuple_features>` or
       :ref:`tuple pair <tuple_pair_features>` features."""

    def __init__(self, mlib):
        BondClasses.__init__(self, mlib, 0)

    def read(self, filename):
        """Read atom tuple information from `filename`.
           This is a text file with a format similar to that accepted by
           :meth:`BondClasses.read`. The file can consist either of sets
           of atom triplets (named with TRPGRP lines and containing triples
           of atoms named on TRIPLET lines) or sets of atom doublets
           using DBLGRP and DOUBLET lines. Each atom but the first in each
           doublet or triplet can also be restricted to match only in
           certain residue types by naming the residue in parentheses before
           the rest of the atom name (and CHARMM-style + or - qualifier).
           For example, a suitable atom triplet file looks like::

             TRPGRP 't1'
               TRIPLET 'ALA' 'CA' '+C' '-C'
             TRPGRP 't2'
               TRIPLET 'ALA' 'CA' '(CYS)+C' '-C'

           The first triplet is named 't1' and will match any set of three
           atoms where the first is called CA in an ALA residue, and the
           other two atoms are C atoms in the previous and next residue.
           The second triplet is similar but will only include triplets where
           the next residue is a CYS.
        """
        return _mdt.mdt_tuple_read(filename, self._mlib._modpt)


class HydrogenBondClasses(BondClasses):
    """Classifications of atoms into hydrogen bond classes.
       Usually accessed as :attr:`Library.hbond_classes`.
       These classes are used by the :class:`features.HydrogenBondAcceptor`,
       :class:`features.HydrogenBondDonor` and
       :class:`features.HydrogenBondSatisfaction` features."""
    def __init__(self, mlib):
        BondClasses.__init__(self, mlib, 1)

    def read(self, filename):
        """Read hydrogen bond atom class information from a file"""
        return _mdt.mdt_hbond_read(filename, self._mlib._modpt)


class TableSection(ModObject):
    """A section of a multi-dimensional table. You should not create
       TableSection objects directly, but rather by indexing a :class:`Table`
       object, as a TableSection is just a 'view' into an existing table.
       For example, ::

         >>> m = mdt.Table(mlib, features=(residue_type, xray_resolution))
         >>> print m[0].entropy()

       would create a section (using m[0]) which is a 1D table over the 2nd
       feature (X-ray resolution) for the first bin (0) of the first feature
       (residue type), and then get the entropy using the
       :meth:`TableSection.entropy` method."""
    _indices = ()
    __mdt = None
    _mlib = None
    _modpt = None
    _basept = None

    def __init__(self, mdt, indices):
        self.__mdt = mdt  # Keep a reference to the MDT
        self._modpt = mdt._modpt
        self._basept = mdt._basept
        self._mlib = mdt._mlib
        self._indices = indices

    def _get_removed_rank(self):
        """Return the number of dimensions removed from the full MDT in this
           section (0 for the full MDT, up to nfeatures - 1)"""
        return len(self._indices)

    def sum(self):
        """Sum of all points in the table"""
        return _mdt.mdt_section_sum(self._basept, self._indices)

    def entropy(self):
        """Entropy of all points in the table"""
        return _mdt.mdt_section_entropy(self._basept, self._indices)

    def mean_stdev(self):
        """Mean and standard deviation of the table"""
        return _mdt.mdt_section_meanstdev(self._basept, self._mlib._modpt,
                                          self._indices)

    def __check_indices(self, indices):
        """Make sure the indices for an MDT section are reasonable"""
        for (feat, indx) in zip(self.features, indices):
            istart, iend = feat.offset, len(feat.bins) + feat.offset - 1
            if indx < 0:
                indx += iend + 1
            if not istart <= indx <= iend:
                raise IndexError("index (%d) not in range %d<=index<=%d"
                                 % (indx, istart, iend))

    def __getitem__(self, indx):
        if isinstance(indx, list):
            indx = tuple(indx)
        elif not isinstance(indx, tuple):
            indx = (indx,)
        if len(indx) < len(self.features):
            self.__check_indices(indx)
            return TableSection(self, self._indices + indx)
        else:
            return _mdt.mdt_get(self._basept, self._indices + indx)

    def __setitem__(self, indx, val):
        if isinstance(indx, list):
            indx = tuple(indx)
        elif not isinstance(indx, tuple):
            indx = (indx,)
        if len(indx) < len(self.features):
            raise ValueError("Cannot set sections of MDTs")
        else:
            _mdt.mdt_set(self._basept, self._indices + indx, val)

    def __get_features(self):
        return _FeatureList(self)

    def __get_offset(self):
        return tuple([f.offset for f in self.features])

    def __get_shape(self):
        return tuple([len(f.bins) for f in self.features])
    features = property(__get_features,
                        doc="Features in this MDT; a list of "
                            ":class:`Feature` objects")
    offset = property(__get_offset,
                      doc="Array offsets; see :attr:`Feature.offset`")
    shape = property(__get_shape, doc="Array shape; the number of "
                                      "bins for each feature")


class Table(TableSection):
    """A multi-dimensional table.

       :Parameters:
         - `mlib`: the MDT `Library` object to use
         - `file`: if specified, the filename to read the initial table from
           (if the name ends with '.hdf5', :meth:`Table.read_hdf5` is used,
           otherwise :meth:`Table.read`)
         - `features`: if specified (and `file` is not), a list of feature
           types to initialize the table with (using :meth:`Table.make`)
         - `bin_type`: type of storage for bin data (see :ref:`binstorage`).
         - `shape`: if specified with `features`, the shape of the new table
           (see :meth:`Table.make`)

       Individual elements from the table can be accessed in standard Python
       fashion, e.g. ::

         >>> import mdt.features
         >>> import modeller
         >>> env = modeller.environ()
         >>> mlib = mdt.Library(env)
         >>> restyp1 = mdt.features.ResidueType(mlib, protein=0)
         >>> restyp2 = mdt.features.ResidueType(mlib, protein=1)
         >>> gap = mdt.features.GapDistance(mlib, mdt.uniform_bins(10, 0, 1))
         >>> m = mdt.Table(mlib, features=(restyp1,restyp2,gap))
         >>> print m[0,0,0]

       You can also access an element as m[0][0][0], a 1D section as m[0][0],
       or a 2D section as m[0]. See :class:`TableSection`.
    """
    _modpt = None
    _basept = None
    _mlib = None
    __views = None

    def __init__(self, mlib, file=None, features=None, bin_type=Double,
                 shape=[]):
        if not isinstance(bin_type, _BinType):
            raise TypeError("bin_type must be a BinType object - "
                            "e.g. mdt.Float, mdt.Double")
        self.__views = []
        self._modpt = _mdt.mdt_new(bin_type._bin_type)
        self._mlib = mlib
        if file:
            if file.endswith(".hdf5"):
                self.read_hdf5(file)
            else:
                self.read(file)
        elif features:
            self.make(features, shape)

    def __getstate__(self):
        d = Table.__getstate__(self)
        d['bin_type'] = _mdt.mod_mdt_bin_type_get(self._basept)
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._modpt = _mdt.mdt_new(d.pop('bin_type'))

    def __del__(self):
        if self._modpt is not None:
            _mdt.mdt_free(self._modpt)

    def __iadd__(self, other):
        _mdt.mdt_add(self._modpt, other._modpt)
        return self

    def __add__(self, other):
        mdtout = self.copy()
        mdtout += other
        return mdtout

    def read(self, file):
        """Read an MDT from `file`.
           ValueError is raised if any views of the table exist
           (see :meth:`Table.get_array_view`)."""
        self.__check_for_views()
        _mdt.mdt_read(self._modpt, self._mlib._modpt, file)

    def read_hdf5(self, file):
        """Read an MDT in HDF5 format from `file`.
           ValueError is raised if any views of the table exist
           (see :meth:`Table.get_array_view`)."""
        self.__check_for_views()
        _mdt.mdt_read_hdf5(self._modpt, self._mlib._modpt, file)

    def __check_for_views(self):
        # Remove any dead weakrefs
        self.__views = [x for x in self.__views if x() is not None]
        if len(self.__views) > 0:
            raise ValueError("Cannot modify the table: views of it exist")

    def __track_view(self, v):
        import weakref
        self.__views.append(weakref.ref(v))

    def get_array_view(self):
        """Get a NumPy array 'view' of this Table. The array contains all of
           the raw data in the MDT table, allowing it to be manipulated with
           NumPy functions. The data are not copied; modifications made to
           the data by NumPy affect the data in the Table (and vice versa).

           Functions that destroy the data in the Table (:meth:`Table.make`,
           :meth:`Table.read` and :meth:`Table.read_hdf5`) cannot be called
           if any NumPy array views exist, since they would invalidate the
           views. The views must first be deleted.

           If MDT was not built with NumPy support, a NotImplementedError
           exception is raised. If NumPy cannot be loaded, an ImportError
           is raised.

           :return: a view of this table.
           :rtype: NumPy array
           """
        v = _mdt.get_numpy(self._modpt, self)
        self.__track_view(v)
        return v

    def copy(self, bin_type=None):
        """
        If `bin_type` is specified, it is the storage type to convert the
        bin data to (see :ref:`binstorage`).

        :return: a copy of this MDT table.
        :rtype: :class:`Table`
        """
        if bin_type is None:
            bin_type = _mdt.mod_mdt_bin_type_get(self._basept)
        elif isinstance(bin_type, _BinType):
            bin_type = bin_type._bin_type
        else:
            raise TypeError("bin_type must be a BinType object - "
                            "e.g. mdt.Float, mdt.Double")
        mdtout = Table(self._mlib)
        _mdt.mdt_copy(self._modpt, mdtout._modpt, bin_type)
        return mdtout

    def make(self, features, shape=[]):
        """Clear the table, and set the features. `features` must be a list of
           previously created objects from the :mod:`mdt.features` module.
           If given, `shape` has the same meaning as in :meth:`Table.reshape`
           and causes the table to use only a subset of the feature bins.

           ValueError is raised if any views of the table exist
           (see :meth:`Table.get_array_view`).
        """
        self.__check_for_views()
        features = self._features_to_ifeat(features)
        _mdt.mdt_make(self._modpt, self._mlib._modpt, features, shape)

    def clear(self):
        """Clear the table (set all bins to zero)"""
        _mdt.mdt_clear(self._modpt)

    def write(self, file, write_preamble=True):
        """Write the table to `file`. If `write_preamble` is False, it will
           only write out the contents of the MDT table, without the preamble
           including the feature list, bins, etc. This is useful for example
           for creating a file to be read by another program, such as
           Mathematica."""
        _mdt.mdt_write(self._modpt, self._mlib._modpt, file, write_preamble)

    def write_hdf5(self, file, gzip=False, chunk_size=1024*1024*10):
        """
        Write an MDT in HDF5 format to `file`.
        Certain library information (such as the mapping from feature
        values to bin indices, and atom or tuple class information)
        and information about the last scan is also written to the file.
        (This information will be missing or incomplete if
        :meth:`add_alignment` hasn't first been called.)
        Note that this information is not read back in by :meth:`read_hdf5`;
        it is intended primarily for other programs that want to reproduce
        the environment in which the MDT was generated as closely as possible.

        :Parameters:
          - `gzip`: If True, compress the table in the HDF5 file with gzip
            using the default compression level; if a number from 0-9, compress
            using that gzip compression level (0=no compression, 9=most);
            if False (the default) do not compress.
          - `chunk_size`: when using gzip, the table must be split up into
            chunks (otherwise it is written contiguously). This parameter
            can either be a list (the same length as the number of features)
            defining the size of each chunk, or it can be the approximate
            number of data points in each chunk, in which case the dimensions
            of the chunk are chosen automatically.
        """
        if gzip is False:
            gzip = -1
        elif gzip is True:
            gzip = 6
        if gzip >= 0 and not hasattr(chunk_size, '__len__'):
            chunk_size = self._guess_chunk_size(self.shape, chunk_size)
        _mdt.mdt_write_hdf5(self._modpt, self._mlib._modpt, file, gzip,
                            chunk_size)

    def _guess_chunk_size(self, shape, chunk_size):
        """Determine a suitable chunk size given the total number of points"""
        import math
        total_size = 1
        for d in shape:
            total_size *= d
        div = math.pow(float(total_size) / float(chunk_size),
                       1. / float(len(shape)))
        ret = []
        for d in shape:
            cd = int(float(d) / div)
            if cd < 1:
                cd = 1
            elif cd > d:
                cd = d
            ret.append(cd)
        return ret

    def reshape(self, features, offset, shape):
        """
        Reorder the MDT features and optionally decrease their ranges.
        When an MDT is created, each feature has exactly the bins defined in
        the `Library`'s bin file. However, for each feature, you can change
        the offset (initial number of bins from the bin file to omit) from the
        default 0, and the shape (total number of bins).

        All parameters should be lists with the same number of elements as
        the MDT has features.

        :Parameters:
          - `features`: the new ordering of the MDT features.
          - `offset`: the new offset (see `offset`).
          - `shape`: the new shape (see `shape`). If any element in this list
            is 0 or negative, it is added to the MDT's existing shape to get
            the new value. Thus, a value of 0 would leave the shape unchanged,
            -1 would remove the last (undefined) bin, etc.
        :return: the reshaped MDT.
        :rtype: :class:`Table`
        """
        features = self._features_to_ifeat(features)
        mdtout = Table(self._mlib)
        _mdt.mdt_reshape(self._modpt, mdtout._modpt, features, offset, shape)
        return mdtout

    def smooth(self, dimensions, weight):
        r"""
        Smooth the MDT with a uniform prior. The MDT is treated either as a
        histogram (if `dimensions` = 1) or a 2D density (`dimensions` = 2)
        of dependent features (the last 1 or 2 features in the table)
        and a uniform distribution is added followed by scaling:

        p\ :sub:`i` = |w1| / n + |w2| |vi| / S

        S = Σ\ :sub:`i`\ :sup:`n` |vi|

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
        :rtype: :class:`Table`

        .. |w1| replace:: w\ :sub:`1`
        .. |w2| replace:: w\ :sub:`2`
        .. |vi| replace:: v\ :sub:`i`
        """
        mdtout = Table(self._mlib)
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
            *p(x | a, b, c, ...)* if `dimensions` = 1, or
            *p(x, y | a, b, c, ...)* if `dimensions` = 2, where y and x are
            the second to last and last features in the list of features.
          - `dx_dy`: widths of the bins (either one or two numbers,
            depending on `dimensions`). If the value of either dx or dy
            is -999, the corresponding bin width is extracted from the MDT
            data structure (not available for all features).
          - `to_zero`: if the histogram is empty, setting this True will set
            the bin values to zero, and False will yield a uniform
            distribution. It has no effect when the histogram is not empty.
          - `to_pdf`: if False, the output is obtained by scaling the input
            such that for 1D histograms Σ :sub:`i` p(x :sub:`i`) = 1,
            and for 2D histograms Σ :sub:`i,j` p(x :sub:`i,j`) = 1. Note
            that `dx_dy` is **not** taken into account during this scaling.

            If it is True, the normalization takes into account `dx_dy` so
            that the normalized distribution is actually a PDF. That is,
            Σ :sub:`i` p(x :sub:`i`) dx = 1 for 1D and
            Σ :sub:`i,j` p(x :sub:`i,j`) dx dy = 1 for 2D, where dx and
            dy are the widths of the bins.
        :return: the normalized MDT.
        :rtype: :class:`Table`

        """
        mdtout = Table(self._mlib)
        _mdt.mdt_normalize(self._modpt, mdtout._modpt, self._mlib._modpt,
                           dimensions, dx_dy, to_zero, to_pdf)
        return mdtout

    def integrate(self, features):
        """
        Integrate the MDT, and reorder the features. This is useful for
        squeezing large MDT arrays into smaller ones, and also for
        eliminating unwanted features (such as X-ray resolution) in
        preparation for :meth:`Table.write`.

        :Parameters:
          - `features`: the new features (all must be present in the
            original MDT).
        :return: the integrated MDT.
        :rtype: :class:`Table`
        """
        features = self._features_to_ifeat(features)
        mdtout = Table(self._mlib)
        _mdt.mdt_integrate(self._modpt, mdtout._modpt, features)
        return mdtout

    def exp_transform(self, offset, expoffset, multiplier, power):
        r"""
        Apply an exponential transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = offset + exp(expoffset + multiplier \* a ^ power)*.

        :rtype: :class:`Table`
        """
        mdtout = self.copy()
        _mdt.mdt_exp_transform(mdtout._basept, offset, expoffset, multiplier,
                               power)
        return mdtout

    def log_transform(self, offset, multiplier, undefined=0.):
        r"""
        Apply a log transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = ln(offset + multiplier \* a)*. Where this would involve the
        logarithm of a negative number, *b* is assigned to be `undefined`.

        :return: the transformed MDT.
        :rtype: :class:`Table`
        """
        mdtout = self.copy()
        _mdt.mdt_log_transform(mdtout._basept, offset, multiplier, undefined)
        return mdtout

    def linear_transform(self, offset, multiplier):
        r"""
        Apply a linear transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = offset + a \* multiplier*.

        :return: the transformed MDT.
        :rtype: :class:`Table`
        """
        mdtout = self.copy()
        _mdt.mdt_linear_transform(mdtout._basept, offset, multiplier)
        return mdtout

    def inverse_transform(self, offset, multiplier, undefined=0.):
        """
        Apply an inverse transform to the MDT.
        Each element in the new MDT, *b*, is obtained from the original
        MDT element *a*, using the following relation:
        *b = offset + multiplier / a*. Where *a* is zero, *b* is
        assigned to be `undefined`.

        :return: the transformed MDT.
        :rtype: :class:`Table`
        """
        mdtout = self.copy()
        _mdt.mdt_inverse_transform(mdtout._basept, offset, multiplier,
                                   undefined)
        return mdtout

    def offset_min(self, dimensions):
        """
        Offset the MDT by the minimum value, either in each 1D section
        (`dimensions` = 1) or in each 2D section (`dimensions` = 2).

        :return: the transformed MDT.
        :rtype: :class:`Table`
        """
        mdtout = self.copy()
        _mdt.mdt_offset_min(mdtout._basept, dimensions)
        return mdtout

    def close(self, dimensions):
        """
        Attempt to 'close' the MDT, so that it is useful for creating splines
        of periodic features.

        If `dimensions` = 1, it makes the two terminal points equal to their
        average. If `dimensions` = 2, it applies the averages to both pairs
        of edges and then again to all four corner points.

        :return: the closed MDT.
        :rtype: :class:`Table`
        """
        mdtout = self.copy()
        _mdt.mdt_close(mdtout._basept, dimensions)
        return mdtout

    def entropy_full(self):
        """Print full entropy information."""
        return _mdt.mdt_entropy_full(self._basept, self._mlib._modpt)

    def entropy_hx(self):
        r"""
        The MDT is integrated to get a 1D histogram, then normalized by
        the sum of the bin values. Finally, entropy is calculated as
        Σ\ :sub:`i` -p\ :sub:`i` ln p\ :sub:`i`

        :return: the entropy of the last dependent variable.
        :rtype: float
        """
        return _mdt.mdt_entropy_hx(self._basept)

    def super_smooth(self, dimensions, prior_weight, entropy_weighing):
        """
        Multi-level smoothing. This super-smoothes the raw frequencies in
        the MDT using the hierarchical smoothing procedure for 1D histograms
        described in Sali and Blundell, JMB 1993. It was also employed in
        Sali and Overington, Prot Sci. 1994.

        Briefly, the idea is to recursively construct the best possible
        prior distribution for smoothing 1D data *p(x | a, b, c, ...)*.
        The best prior is a weighted sum (weights optionally based on
        entropy) of the best possible estimate of *p(x | a, b, ...)*
        integrated over *c* for each *c*. Each one of these can itself be
        obtained from a prior and the data, and so on recursively.

        The example above is for a single dependent feature (*x*), which is the
        case when `dimensions` = 1. *x* should be the last feature in the
        table. `dimensions` can be set to other values if you have more
        dependent features - for example, `dimensions` = 2 will work with
        *p(x, y | a, b, c, ...)* where *x* and *y* are the last two features
        in the table.

        :Parameters:
          - `dimensions`: Number of dependent features.
          - `prior_weight`: Weight for the prior distribution.
          - `entropy_weighing`: Whether to weight distributions by their
            entropies.
        :return: the smoothed MDT.
        :rtype: :class:`Table`
        """
        mdtout = Table(self._mlib)
        _mdt.mdt_super_smooth(self._modpt, mdtout._modpt, dimensions,
                              prior_weight, entropy_weighing)
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
        return _mdt.mdt_write_asgl(self._basept, self._mlib._modpt, asglroot,
                                   text, dimensions, every_x_numbered,
                                   every_y_numbered, plot_density_cutoff,
                                   plots_per_page, plot_position, plot_type,
                                   x_decimal, y_decimal)

    def add_alignment(self, aln, distngh=6.0, surftyp=1, accessibility_type=8,
                      residue_span_range=(-99999, -2, 2, 99999),
                      chain_span_range=(-99999, 0, 0, 99999),
                      bond_span_range=None, disulfide=False,
                      exclude_bonds=False, exclude_angles=False,
                      exclude_dihedrals=False, sympairs=False,
                      symtriples=False, io=None, edat=None):
        """
        Add data from a Modeller alignment to this MDT.
        This method will first scan through all proteins, pairs of proteins,
        or triples of proteins in the alignment (it will scan all triples if
        the :class:`mdt.Library` contains features defined on all of
        proteins 0, 1 and 2, pairs if the features are defined on two
        different proteins, and individual proteins otherwise). Within each
        protein, it may then scan through all residues, atoms, etc. if the
        features request it (see :ref:`the scan types table <scantypes>`).

        :Parameters:
          - `aln`: Modeller alignment.
          - `distngh`: distance below which residues are considered neighbors.
            Used by :class:`features.NeighborhoodDifference`.
          - `surftyp`: 1 for PSA contact area, 2 for surface area.
            Used by :class:`features.AtomAccessibility`.
          - `accessibility_type`: PSA accessibility type (1-10).
            Used by :class:`features.AtomAccessibility`.
          - `residue_span_range`: sequence separation (inclusive) for
            :ref:`residue pair <residue_pair_features>`,
            :ref:`atom pair <atom_pair_features>` and
            :ref:`tuple pair <tuple_pair_features>` features. For the two
            residue indices r1 and r2 in the tuple-tuple and atom- atom cases,
            or two alignment position indices in the residue-residue case,
            the following must be true:

            *residue_span_range[0] <= (r2 - r1) <= residue_span_range[1]*

            *residue_span_range[2] <= (r2 - r1) <= residue_span_range[3]*

            For symmetric residue-residue features, only one condition
            must be met:

            *residue_span_range[2] <= abs(r2 - r1) <= residue_span_range[3]*

            For example, the default value of (-99999, -2, 2, 99999) excludes
            all pairs within the same residue (for which the sequence
            separation is 0) or within adjacent residues (for which the
            separation is 1 or -1).
          - `chain_span_range`: works like `residue_span_range`, but for the
            chain indices. It is used only by the
            :ref:`atom pair <atom_pair_features>` and
            :ref:`tuple pair <tuple_pair_features>` features. The default value
            of (-99999, 0, 0, 99999) allows all interactions. For example,
            using (-99999, -1, 1, 99999) instead would exclude all interactions
            within the same chain.
          - `bond_span_range`: if given, it should be a list of two integers
            which specify the minimum and maximum number of bonds that separate
            a pair of atoms in the scan. It is used only by the
            :ref:`atom pair <atom_pair_features>` and
            :ref:`tuple pair <tuple_pair_features>` features. (See
            :class:`features.AtomBondSeparation` for more details.) The bond
            library (see :attr:`Library.bond_classes`) must be loaded to use
            this. For example, using (1, 2) will include only atoms that
            are directly chemically bonded or that are both bonded to a third
            atom, while (0, 9999) will only exclude pairs of atoms that have
            no path of bonds between them (e.g. atoms in different chains or
            when at least one of the atoms is not involved in any bonds).
            As a special case, if the maximum span is negative, no limit is
            enforced. For example, (2, 99999) will include all atoms that have
            a path of bonds between them except directly bonded pairs (and
            thus exclude pairs in different chains) while (2, -1) will also
            include inter-chain interactions.
          - `disulfide`: if True, then the `bond_span_range` considers
            disulfide bonds (defined as any pair of SG atoms in CYS residues
            less than 2.5 angstroms apart) when calculating the bond separation
            between atoms. Only disulfide bridges within 3 residues of the
            atom pair are considered for computational efficiency.
          - `exclude_bonds`: if True, then all pairs of atoms involved in a
            chemical bond (see :attr:`Library.bond_classes`) are excluded from
            :ref:`atom pair <atom_pair_features>` and
            :ref:`tuple pair <tuple_pair_features>` features.
          - `exclude_angles`: if True, then the 1-3 pair of atoms from each
            angle are excluded (see `exclude_bonds`).
          - `exclude_dihedrals`: if True, then the 1-4 pair of atoms from each
            dihedral are excluded (see `exclude_bonds`).
          - `sympairs`: if True, then protein pair scans are done in a
            symmetric fashion - e.g. when scanning an alignment of A, B and
            C, the following pairs are scanned: AB, BC, AC. By default a
            non-symmetric scan is performed, scanning AB, BC, AC, BA, CB, CA.
          - `symtriples`: if True, then protein triple scans are done in a
            symmetric fashion - e.g. when scanning an alignment of A, B and
            C, the following triples are scanned: ABC, ACB, BAC. By default a
            non-symmetric scan is performed, scanning ABC, ACB, BAC, CBA,
            BCA, CAB.
        """
        if io is None:
            io = self._mlib._env.io
        if edat is None:
            edat = self._mlib._env.edat
        _mdt.mdt_add_alignment(self._modpt, self._mlib._modpt, aln.modpt,
                               distngh, False, surftyp, accessibility_type,
                               residue_span_range, chain_span_range,
                               _prepare_bond_span(bond_span_range), disulfide,
                               exclude_bonds, exclude_angles,
                               exclude_dihedrals, sympairs, symtriples,
                               io.modpt, edat.modpt)

    def add_alignment_witherr(self, aln, distngh=6.0, surftyp=1,
                              accessibility_type=8,
                              residue_span_range=(-99999, -2, 2, 99999),
                              chain_span_range=(-99999, 0, 0, 99999),
                              bond_span_range=None, disulfide=False,
                              exclude_bonds=False, exclude_angles=False,
                              exclude_dihedrals=False,
                              sympairs=False, symtriples=False, io=None,
                              edat=None, errorscale=1):
        """
        Add data from a Modeller alignment to this MDT. Same as add_alignment
        except the errors in data are taken into account.
        The parameter errorscale controls how the error is used:

          - `0`: the errors are ignored; this function is the same as
                 add_alignment.
          - `>0` : the errors are taken into account by propagating the errors
                  in each axis of each atom into the calculated distances
                  or angles. The errors in the position of individual
                  atoms are first calculated using B-iso, X-ray resolution,
                  and R-factor, and then divided by this errorscale value.
        """
        if io is None:
            io = self._mlib._env.io
        if edat is None:
            edat = self._mlib._env.edat
        _mdt.mdt_add_alignment_witherr(self._modpt, self._mlib._modpt,
                                       aln.modpt, distngh, False, surftyp,
                                       accessibility_type, residue_span_range,
                                       chain_span_range,
                                       _prepare_bond_span(bond_span_range),
                                       disulfide,
                                       exclude_bonds,
                                       exclude_angles, exclude_dihedrals,
                                       sympairs, symtriples, io.modpt,
                                       edat.modpt, errorscale)

    def open_alignment(self, aln, distngh=6.0, surftyp=1, accessibility_type=8,
                       sympairs=False, symtriples=False, io=None, edat=None):
        """
        Open a Modeller alignment to allow MDT indices to be queried
        (see :class:`Source`). Arguments are as for
        :meth:`Table.add_alignment`.

        :rtype: :class:`Source`
        """
        return Source(self, self._mlib, aln, distngh, surftyp,
                      accessibility_type, sympairs, symtriples, io, edat)

    def _features_to_ifeat(self, features):
        """Utility function to map objects from `mdt.features` to
           integer feature types"""
        ifeat = []
        if not isinstance(features, (list, tuple)):
            features = (features,)
        for feat in features:
            if hasattr(feat, '_get_ifeat'):
                ifeat.append(feat._get_ifeat(self._mlib))
            else:
                raise TypeError("features should be objects from mdt.features")
        return ifeat

    def __get_pdf(self):
        return _mdt.mdt_pdf_get(self._modpt)

    def __get_n_proteins(self):
        return _mdt.mdt_n_proteins_get(self._modpt)

    def __get_n_protein_pairs(self):
        return _mdt.mdt_n_protein_pairs_get(self._modpt)

    def __get_sample_size(self):
        return _mdt.mdt_sample_size_get(self._modpt)

    def __get_symmetric(self):
        return _mdt.mdt_symmetric_get(self._modpt)

    def __get_basept(self):
        return _mdt.mdt_base_get(self._modpt)

    _basept = property(__get_basept)
    pdf = property(__get_pdf, doc="Whether this MDT is a PDF")
    n_proteins = property(__get_n_proteins, doc="Number of proteins")
    n_protein_pairs = property(__get_n_protein_pairs,
                               doc="Number of protein pairs")
    sample_size = property(__get_sample_size, doc="Number of sample points")
    symmetric = property(__get_symmetric,
                         doc="True if a symmetric scan can be performed")


class _FeatureList(modlist.FixList):
    """A list of all features in an MDT."""

    def __init__(self, mdt):
        self.__mdt = mdt
        self.__removed_rank = mdt._get_removed_rank()
        modlist.FixList.__init__(self)

    def __len__(self):
        return _mdt.mod_mdt_nfeat_get(self.__mdt._basept) - self.__removed_rank

    def _getfunc(self, indx):
        return Feature(self.__mdt, indx + self.__removed_rank)


class Feature(object):
    """A single feature in an MDT. Generally accessed as
       :attr:`TableSection.features`."""

    def __init__(self, mdt, indx):
        self._mdt = mdt
        self._indx = indx

    def __get_ifeat(self):
        return _mdt.mod_mdt_feature_ifeat_get(self._modpt)

    def __get_bins(self):
        return _BinList(self)

    def __get_offset(self):
        return _mdt.mod_mdt_feature_istart_get(self._modpt) - 1

    def __get_periodic(self):
        return _mdt.mdt_feature_periodic_get(self._mdt._mlib._modpt,
                                             self.ifeat)

    def __get_modpt(self):
        return _mdt.mod_mdt_feature_get(self._mdt._basept, self._indx)

    _modpt = property(__get_modpt)
    ifeat = property(__get_ifeat, doc="Integer type")
    bins = property(__get_bins,
                    doc="Feature bins; a list of :class:`Bin` objects")
    offset = property(__get_offset,
                      doc="Offset of first bin compared to the MDT library "
                          "feature (usually 0, but can be changed with "
                          ":meth:`Table.reshape`)")
    periodic = property(__get_periodic, doc="Whether feature is periodic")


class _BinList(modlist.FixList):
    """A list of all bins in a feature."""

    def __init__(self, feature):
        self.__feature = feature
        self._mdt = feature._mdt
        modlist.FixList.__init__(self)

    def __len__(self):
        return _mdt.mod_mdt_feature_nbins_get(self.__feature._modpt)

    def _getfunc(self, indx):
        return Bin(self.__feature, indx)


class Bin(object):
    """A single bin in a feature. Generally accessed as
       :attr:`Feature.bins`."""

    def __init__(self, feature, indx):
        self.__feature = feature
        self.__indx = indx

    def __get_symb(self):
        return _mdt.mod_mdt_bin_symbol_get(self._modpt)

    def __get_range(self):
        return (_mdt.mod_mdt_bin_rang1_get(self._modpt),
                _mdt.mod_mdt_bin_rang2_get(self._modpt))

    def __get_modpt(self):
        nfeat = self.__feature._indx
        mdt = self.__feature._mdt._basept
        mlib = self.__feature._mdt._mlib._modpt
        return _mdt.mdt_library_bin_get(mdt, mlib, nfeat, self.__indx)

    _modpt = property(__get_modpt)
    symbol = property(__get_symb, doc="Bin symbol")
    range = property(__get_range,
                     doc="Bin range; usually the minimum and maximum "
                         "floating-point values for the feature to be "
                         "placed in this bin.")


class Source(object):
    """A source of data for an MDT (generally a Modeller alignment, opened
       with :meth:`Table.open_alignment`)."""

    def __init__(self, mdt, mlib, aln, distngh, surftyp, accessibility_type,
                 sympairs, symtriples, io, edat):
        self._mdt = mdt
        self._mlib = mlib
        self._aln = aln
        if io is None:
            io = mlib._env.io
        if edat is None:
            edat = mlib._env.edat
        self._edat = edat
        self._modpt = _mdt.mdt_alignment_open(mdt._modpt, mlib._modpt,
                                              aln.modpt, distngh, False,
                                              surftyp, accessibility_type,
                                              sympairs, symtriples, io.modpt)

    def __del__(self):
        if hasattr(self, "_modpt"):
            _mdt.mdt_alignment_close(self._modpt)

    def sum(self, residue_span_range=(-99999, -2, 2, 99999),
            chain_span_range=(-99999, 0, 0, 99999),
            bond_span_range=None, disulfide=False,
            exclude_bonds=False, exclude_angles=False,
            exclude_dihedrals=False):
        """Scan all data points in the source, and return the sum.
           See :meth:`Table.add_alignment` for a description of the
           `residue_span_range`, `chain_span_range` and `exclude_*`
           arguments."""
        f = _mdt.mdt_source_sum
        return f(self._modpt, self._mdt._modpt, self._mlib._modpt,
                 residue_span_range, chain_span_range,
                 _prepare_bond_span(bond_span_range), disulfide,
                 exclude_bonds, exclude_angles, exclude_dihedrals,
                 self._edat.modpt)

    def index(self, feat, is1, ip1, is2, ir1, ir2, ir1p, ir2p, ia1, ia1p,
              ip2, ibnd1, ibnd1p, is3, ir3, ir3p):
        """
        Return the bin index (starting at 1) of a single MDT feature.
        (Arguments ending in 2 and 3 are used for features involving pairs
        or triples of proteins.)

        .. warning:: This is a low-level interface, and no bounds checking is
           performed on these parameters. Avoid this function if possible.

        :Parameters:
          - `feat`: MDT feature object from :mod:`mdt.features` module.
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
        return f(self._modpt, feat._ifeat, is1, ip1, is2, ir1, ir2, ir1p, ir2p,
                 ia1, ia1p, ip2, ibnd1, ibnd1p, is3, ir3, ir3p,
                 self._mlib._modpt, self._edat.modpt)


def _pass_cutoffs(mdt, num, bin, density_cutoff, entropy_cutoff):
    if density_cutoff is not None:
        sum = mdt[num].sum()
        if sum < density_cutoff:
            print("Restraint %s skipped: density %.4f below cutoff %.4f"
                  % (bin.symbol, sum, density_cutoff))
            return False
    if entropy_cutoff is not None:
        entropy = mdt[num].entropy()
        if entropy > entropy_cutoff:
            print("Restraint %s skipped: entropy %.4f above cutoff %.4f"
                  % (bin.symbol, entropy, entropy_cutoff))
            return False
    return True


def _write_meanstdevlib(fh, mdt, numat, phystype, feattype, convfunc,
                        density_cutoff=None, entropy_cutoff=None):
    fh.write("#   residue    atoms        mean    stdev\n")
    fh.write("_params = [\n")
    for (num, bin) in enumerate(mdt.features[0].bins):
        if _pass_cutoffs(mdt, num, bin, density_cutoff, entropy_cutoff):
            symbols = bin.symbol.split(':')
            res = symbols[0]
            ats = tuple(symbols[1:])
            if len(ats) != numat:
                example_atoms = ['CA', 'C', 'N', 'O']
                raise ValueError("Bin name '%s' should be residue name "
                                 "plus %d atoms, separated by colons, e.g. %s"
                                 % (bin.symbol, numat,
                                    'ALA:' + ':'.join(example_atoms[:numat])))
            mean, stdev = mdt[num].mean_stdev()
            fh.write("    ( '%s', %s, %.4f, %.4f ),\n"
                     % (res, str(ats), convfunc(mean), convfunc(stdev)))
    fh.write("""  ]

def make_restraints(atmsel, restraints, num_selected):
    from modeller import forms, physical, features
    for (res, atoms, mean, stdev) in _params:
        for a in atmsel.find_atoms(res, atoms, num_selected):
            r = forms.gaussian(%s, %s(*a), mean,
                               stdev)
            restraints.add(r)\n""" % (phystype, feattype))


def _noconv(a):
    return a


def _degrees_to_radians(a):
    pi = 3.14159265358979323846
    return a / 180.0 * pi


def write_bondlib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    """
    Write out a Modeller bond library file from an MDT. The input MDT should be
    a 2D table (usually of bond type and bond distance). For each bond type,
    the 1D MDT section (see :class:`TableSection`) of bond distance is
    examined, and its mean and standard deviation used to generate a
    Modeller harmonic restraint.

    :Parameters:
      - `fh`: Python file to write to
      - `mdt`: input MDT :class:`Table` object
      - `density_cutoff`: if specified, MDT bond distance sections with sums
        below this value are not used
      - `entropy_cutoff`: if specified, MDT bond distance sections with
        entropies above this value are not used
    """
    _write_meanstdevlib(fh, mdt, 2, "physical.bond", "features.distance",
                        _noconv, density_cutoff, entropy_cutoff)


def write_anglelib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    """
    Write out a Modeller angle library file from an MDT. See
    :func:`write_bondlib` for more details. The MDT should be a 2D table,
    usually of angle type and bond angle.
    """
    _write_meanstdevlib(fh, mdt, 3, "physical.angle", "features.angle",
                        _degrees_to_radians, density_cutoff, entropy_cutoff)


def write_improperlib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    """
    Write out a Modeller dihedral angle library file from an MDT. See
    :func:`write_bondlib` for more details. The MDT should be a 2D table,
    usually of dihedral type and bond dihedral angle.
    """
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


def write_statpot(fh, mdt):
    """
    Write out a Modeller statistical potential file (as accepted by
    group_restraints.append()). The MDT is assumed to be a 3D table of distance
    against the types of the two atoms. No special processing
    is done, so it is expected that the user has first done any necessary
    transformations (e.g. normalization with :meth:`Table.normalize` to
    convert raw counts into a PDF, negative log transform with
    :meth:`Table.log_transform` and :meth:`Table.linear_transform` to
    convert a PDF into a statistical potential).
    """
    modality = len(mdt.features[0].bins)
    delta = mdt.features[0].bins[0].range[1]
    low = mdt.features[0].bins[0].range[0] + delta / 2.
    high = mdt.features[0].bins[-1].range[0] + delta / 2.
    fh.write("MOD5\n")
    fmt = "R" + " %4d" * 7 + " %5s   %5s  " + " %9.4f" * 6 + " "
    for na1, a1 in enumerate(mdt.features[1].bins):
        for na2, a2 in enumerate(mdt.features[2].bins[na1:]):
            splinevals = ["%9.4f" % mdt[x, na1, na2] for x in range(modality)]
            fh.write(fmt % (10, modality, 1, 31, 2, modality + 6, 0, a1.symbol,
                            a2.symbol, 0., low, high, delta, 0., 0.)
                     + " ".join(splinevals) + "\n")


def write_splinelib(fh, mdt, dihtype, density_cutoff=None,
                    entropy_cutoff=None):
    """
    Write out a Modeller 1D spline library file from an MDT.
    The MDT should be a 2D table, usually of residue type and a chi dihedral
    angle. `dihtype` should identify the dihedral type
    (i.e. chi1/chi2/chi3/chi4). The operation is similar to
    :func:`write_bondlib`,
    but each MDT section is treated as the spline values. No special processing
    is done, so it is expected that the user has first done any necessary
    transformations (e.g. normalization with :meth:`Table.normalize` to
    convert raw counts into a PDF, negative log transform with
    :meth:`Table.log_transform` and :meth:`Table.linear_transform` to
    convert a PDF into a statistical potential).
    """
    (periodic, dx, x1, x2) = _get_splinerange(mdt.features[1])

    fh.write("#   residue   spline values\n")
    fh.write("_params = [\n")
    for (nx, bin) in enumerate(mdt.features[0].bins):
        if _pass_cutoffs(mdt, nx, bin, density_cutoff, entropy_cutoff):
            splinevals = ["%.4f" % mdt[nx, ny]
                          for ny in range(len(mdt.features[1].bins))]
            fh.write("    ( '%s', (%s) ),\n"
                     % (bin.symbol, ', '.join(splinevals)))
    fh.write("""  ]

def make_restraints(atmsel, restraints, num_selected):
    from modeller import forms, physical, features
    for (res, values) in _params:
        arr = True
        for a in atmsel.find_%s_dihedrals(res, num_selected):
            r = forms.spline(physical.%s_dihedral,
                             features.dihedral(*a), open=%s, low=%.5f,
                             high=%.5f, delta=%.5f, lowderiv=0,
                             highderiv=0, values=values, use_array=arr)
            arr = restraints.add(r)\n""" % (dihtype, dihtype, not periodic,
                                            x1, x2, dx))


def write_2dsplinelib(fh, mdt, density_cutoff=None, entropy_cutoff=None):
    """
    Write out a Modeller 2D spline library file from an MDT.
    See :func:`write_splinelib` for more details. The input MDT should be
    a 3D table, e.g. of residue type, phi angle, and psi angle.
    """
    (yperiodic, dy, y1, y2) = _get_splinerange(mdt.features[1])
    (zperiodic, dz, z1, z2) = _get_splinerange(mdt.features[2])

    fh.write("#   residue   spline values\n")
    fh.write("_params = [\n")
    for (nx, bin) in enumerate(mdt.features[0].bins):
        if _pass_cutoffs(mdt, nx, bin, density_cutoff, entropy_cutoff):
            splinevals = []
            for ny in range(len(mdt.features[1].bins)):
                for nz in range(len(mdt.features[2].bins)):
                    splinevals.append("%.4f" % mdt[nx, ny, nz])
            fh.write("    ( '%s', (%s) ),\n"
                     % (bin.symbol, ', '.join(splinevals)))
    fh.write("""  ]

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
            arr = restraints.add(r)\n""" % (len(mdt.features[1].bins),
                                            len(mdt.features[2].bins),
                                            not yperiodic, y1, y2, dy,
                                            not zperiodic, z1, z2, dz))


def uniform_bins(num, start, width):
    """Make a list of `num` equally-sized bins, each of which has the given
       `width`, and starting at `start`. This is suitable for input to any of
       the classes in :mod:`mdt.features` which need a list of bins."""
    bins = []
    for i in range(num):
        st = start + width * i
        en = st + width
        sym = "%.1f" % st
        bins.append((st, en, sym))
    return bins
