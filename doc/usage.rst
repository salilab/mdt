.. highlight:: rest

.. currentmodule:: mdt

Usage
=====

.. _features:

MDT features
------------

A 'feature' in MDT is simply some binnable property of your input alignment.
Example features include the
:class:`residue type <features.ResidueType>`,
:class:`chi1 <features.Chi1Dihedral>`
and :class:`Phi <features.PhiDihedral>` dihedral angles,
:class:`sequence identity <features.SequenceIdentity>` between two sequences,
:class:`X-ray resolution <features.XRayResolution>`,
:class:`atom-atom distances <features.AtomDistance>`,
:class:`atom type <features.AtomType>`, and
:class:`bond length <features.BondLength>`.

MDT understands that different features act on different sets of proteins,
or parts of proteins, and will automatically scan over the correct range to
collect necessary statistics (e.g. when you call :meth:`Table.add_alignment`).
For example, to collect statistics for the residue type feature, it is
necessary to scan all residues in all proteins in the alignment. The
X-ray resolution feature, on the other hand, only requires each protein
in the alignment to be scanned, not each residue.
The atom-atom distance feature requires scanning over all pairs of atoms in
all proteins in the alignment, while the sequence identity feature requires
scanning all pairs of proteins in the alignment. If you construct a table of
multiple features, the most fine-grained of the features determines the scan -
for example, a table of X-ray resolution against &Phi; dihedral would require
a scan of all residues.
See :ref:`the scan types table <scantypes>` for all of the scan types.

When choosing which proteins to scan, MDT also considers the features. It
will scan each protein individually, all pairs of proteins, or all triples of
proteins. The latter two scans only happen if you have features in your table
that require multiple proteins (e.g.
:class:`features.ProteinPair` or :class:`features.AlignedResidue`
features) or you have single-protein features such as
:class:`features.Protein` or :class:`features.Residue`
but you have asked to evaluate them on the second or third protein (by setting
the *protein* argument to 1 or 2 rather than the default 0).

MDT also knows that some
:class:`residue pair <features.ResiduePair>` or
:class:`atom pair <features.AtomPair>` features are symmetric,
and will perform a non-redundant scan in this case. If, however, any feature
in the table is asymmetric, a full scan is performed. If in doubt, you can query
:attr:`Table.symmetric` to see whether
a symmetric scan will be performed for the current set of features.
(Currently, any :class:`tuple pair <features.TuplePair>`
feature in your table forces a full scan.)

The feature bins determine how to convert a feature value into a frequency
table.
For most feature types, you can specify how many bins to use, and their value
ranges - see :ref:`binspec` for more information. The last bin
is always reserved as an 'undefined' bin, for values that don't fall into
any other bin [#undefbin]_.

(Some features are predetermined by the setup of the system - for
example, the :class:`residue type <features.ResidueType>`
feature always has 22 bins - 20 for the standard amino acids, 1 for gaps in
the alignment, and 1 for undefined.)

.. _scantypes:

======================================  ========================================
Type                                    Example feature
======================================  ========================================
Protein                                 :class:`features.XRayResolution`
Residue [#resscan]_                     :class:`features.Chi1Dihedral`
Residue pair [#resscan]_ [#pairscan]_   :class:`features.ResidueIndexDifference`
Atom                                    :class:`features.AtomType`
Atom pair [#pairscan]_                  :class:`features.AtomDistance`
Atom tuple                              :class:`features.TupleType`
Atom tuple pair                         :class:`features.TupleDistance`
Chemical bond                           :class:`features.BondType`
Chemical angle                          :class:`features.Angle`
Chemical dihedral angle                 :class:`features.Dihedral`
======================================  ========================================

.. _depfeat:

Dependent and independent features
----------------------------------

An MDT :class:`Table` object is simply a table of
counts N(a,b,c,...,d) for features a, b, c, ..., d. However, this is often used
to generate a conditional PDF, p(x,y,...,z | a,b,...,c) for independent features
a, b, ..., c and dependent features x, y, ..., z. By convention in MDT the
dependent features are the last or rightmost features in the table, and so
methods which are designed to deal with PDFs such as
:meth:`Table.smooth`, :meth:`Table.super_smooth`,
:meth:`Table.normalize`, :meth:`Table.offset_min`, :meth:`Table.close`
expect the dependent features to be the last features. If necessary you can
reorder the features using :meth:`Table.integrate`.

.. _binspec:

Specification of bins
---------------------

Most features take a *bins* argument when they are created,
which specifies the bin ranges. This is simply a list of (start, end, symbol)
triples, which specify the feature range for each bin, and the symbol to refer
to it by. For example, the following creates an
:class:`X-ray resolution <features.XRayResolution>` feature,
with 4 bins, the first for 0.51-1.4 |Angstrom|,
the second for 1.4-1.6 |Angstrom|, and so on. Anything below 0.51 |Angstrom| or 
2.0 |Angstrom| or above (or an undefined value) will be placed into a fifth
'undefined' bin.

.. code-block:: python

   xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 1.4, "<1.4"),
                                                  (1.4,  1.6, "1.4-1.6"),
                                                  (1.6,  1.8, "1.6-1.8"),
                                                  (1.8,  2.0, "1.8-2.0")])

.. note::
   Bin ranges in MDT are half-closed, i.e. a feature value must be greater
   than or equal to the lower value of the range, and less than the upper
   value, to be counted in the bin. For example, in the case above,
   1.0 |Angstrom| would be placed into the first bin, and 1.4 |Angstrom|
   into the second. (If you define bins with overlapping ranges, values
   will be placed into the first bin that matches.)

In most cases, a set of bins of equal width is desired, and it is
tedious to specify these by hand. A utility function
:func:`uniform_bins` is provided, which takes
three arguments - the number of bins, the lower range of the first bin,
and the width of each bin - and creates a set of bins; all bins are of the
same size and follow after the first bin. For example, the following bins the
:class:`atom-atom distance <features.AtomDistance>` feature into 60 bins,
each 0.5 |Angstrom| wide, with the first bin starting at 0 |Angstrom|.
The first bin is thus 0-0.5 |Angstrom|, the second 0.5-1.0 |Angstrom|, and so
on, up to bin 60 which is 29.5-30.0 |Angstrom|. The additional 'undefined'
bin thus counts anything below 0 |Angstrom|, greater than or equal to
30.0 |Angstrom|, or which could not be calculated for some reason).

.. code-block:: python

   atdist = mdt.features.AtomDistance(mlib, bins=mdt.uniform_bins(60, 0, 0.5))

.. _binstorage:

Storage for bin data
--------------------

By default, when a table is created in MDT it uses double precision floating
point to store the counts. This allows large counts themselves to be accurately
scored, and can also store floating point data such as PDFs. However, for
very large tables, this may use a prohibitive amount of memory. Therefore, it
is possible to change the data type used to store bin data, by specifying
the *bin_type* parameter when creating a
:class:`Table` object. The same parameter can be
given to :meth:`Table.copy`, to make a copy
of the table using a different data type for its storage. Note that other
data types use less storage, but can also store a smaller range of counts.
For example, the :class:`UnsignedInt8`
data type uses only a single byte for each bin, but can only store integer
counts between 0 and 255 (floating point values, or values outside of this
range, will be truncated). MDT uses double precision floating point for all
internal operations, but any storage of bin values uses the user-selected
bin type. Thus you should be careful not to use an inappropriate bin type -
for example, don't use an integer bin type if you are planning to store PDFs
or perform normalization, smoothing, etc.

.. _compilation:

Compilation from source code
============================

You can get the current MDT code by running the following::

   svn co https://svn.salilab.org/impmod/trunk impmod

If you already have a copy of MDT, you can update it to the current code
by running::

   svn update

To compile, create a file called :file:`config.py` in the
:file:`impmod/mdt` directory, and in it set the
*modeller* Python variable to the directory where you have
MODELLER installed. Then run `scons` in the same directory
(and optionally `scons test`). To install, run
`scons prefix=/foo install`, which will install MDT
in the :file:`/foo` directory.

.. _running:

Example MDT script
==================

MDT is simply a Python extension module, and as such can be used in
combination with other Python modules, such as MODELLER.

Generally speaking, to use MDT, you should first create a
:class:`Library` object. You can then read in any
necessary additional files, such as
the definitions of chemical bonds (see :ref:`chembonds` for an
example), or atom tuples. Next, you can define one or more features, which are
classes in the :mod:`features` module. Finally, you can create one or more
:class:`Table` objects, using a selection of the
features you added to the Library, to hold the frequency tables
themselves. See :ref:`chi1` for an example.

To run a script :file:`foo.py`, simply run::
   ${MODINSTALLSVN}/bin/modpy.sh python foo.py


.. rubric:: Footnotes

.. [#undefbin] You can, however, remove the 'undefined' bin using
   :meth:`Table.reshape`.

.. [#resscan] Residue and residue pair scans are also used for
   'one atom per residue' features, such as
   :class:`features.ResidueDistance`, which is the distance between the
   'special atom' in two residues. This special atom is usually
   C&alpha, but can be overridden by specifying the
   *distance_atoms* parameter when creating the
   :class:`Library` object.

.. [#pairscan] When looking at pairs of atoms or residues, it is useful to
   extract information about the 'other' atom or residue in the pair. This
   other atom or residue is termed 'pos2' in MDT, and can be asked for when
   creating the feature. For example, when building a table of atom-atom
   distances (:class:`features.AtomDistance`) it may be useful to tabulate
   it against the atom types of both the first atom. This is done by also
   using two copies of the :class:`~features.AtomType`, the second with
   *pos2=True*.

.. |Angstrom| unicode:: U+212B
