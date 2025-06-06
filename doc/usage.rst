.. highlight:: rest

.. currentmodule:: mdt

Usage
=====

MDT is simply a Python extension module, and as such can be used in
combination with other Python modules, such as MODELLER or the Python
standard library.

Running pre-built binaries
--------------------------

The easiest way to use MDT is to install the pre-built binary RPM for your
variety of Linux (this will first require you to install the Modeller RPM).
Then you should simply be able to run an MDT script :file:`foo.py` just like
any regular Python script with a command similar to::

   python3 foo.py

In the Sali lab, MDT is built as part of the nightly build system at the
same time as MODELLER. Thus you can set up your system to run MDT scripts
by running::

   module load modeller

.. _anaconda:

Using with Anaconda Python
--------------------------

There is an MDT package available for `Anaconda Python <https://store.continuum.io/cshop/anaconda/>`_ for Mac and Linux.
To install it, simply run::

    conda install -c salilab mdt

.. _homebrew:

Using a Homebrew package
------------------------

If you are using a Mac with `homebrew <https://brew.sh>`_ you can
get MDT by running in a terminal window::

    brew tap salilab/salilab; brew install salilab/salilab/mdt

If you don't already have Modeller installed, you can get it by running
`brew install modeller` before you install MDT.

.. _compilation:

Compilation from source code
----------------------------

The MDT source code can be downloaded from
`GitHub <https://github.com/salilab/mdt/>`_.

Install dependent packages needed for MDT: MODELLER, glib, SWIG,
pkg-config, Python, and HDF5:

* MODELLER 10.0 or later is required.
* glib 2.4 or later is required. It is available as pre-built packages for most
  modern Linux distributions; there is also a MacPorts package for Mac users.
* SWIG 1.3.39 or later is required.
* Unfortunately HDF5 only works if you use the exact same version that is
  used by MODELLER. See the MODELLER ChangeLog for the version to use.
* Python 2.7, or any version of Python 3, is required.

To compile, run `scons` in the same directory (and optionally `scons test`)
to build (and test) MDT. This will produce a script :file:`bin/mdtpy.sh`
which can be used to run an MDT Python script :file:`foo.py`::

   bin/mdtpy.sh python3 foo.py

.. note::
   If you didn't use the RPM or Debian package to install Modeller then you
   will need to tell MDT where it can find Modeller. To do this, create a file
   called :file:`config.py`, and in it set the *modeller* Python variable to
   the directory where you have MODELLER installed (on a Mac, this would look
   like `modeller="/Library/modeller-XXX"` where `XXX` is the Modeller version).

   If you installed any of the prerequisites in non-standard locations (i.e.
   not `/usr/include` for glib and HDF5, and not `/usr/bin` for pkg-config
   or SWIG) you will also need to tell scons where to find them. Add similar
   lines to :file:`config.py` to set *path* for pkg-config and SWIG
   and *includepath* for glib and HDF5 (e.g. `path="/opt/local/bin"`
   and `includepath="/opt/local/include"` on a Mac).

If you want to install MDT, run `scons install`. You can additionally specify
a `prefix` option (or set it in :file:`config.py`) to install in a different
directory. For example, `scons prefix=/foo install` will install MDT in
the :file:`/foo` directory.

.. _running:

Example MDT script
------------------

Generally speaking, to use MDT, you should

 #. Create a :class:`Library` object.
 #. Read any necessary additional files into the library, such as the
    definitions of chemical bonds (see :ref:`chembonds` for an example),
    or atom tuples.
 #. Define one or more features, which are classes in the :mod:`mdt.features`
    module.
 #. Create one or more :class:`Table` objects, using a selection of the
    features you added to the Library, to hold the frequency tables
    themselves.
 #. Collect statistics into the table using methods such as
    :meth:`Table.add_alignment`.
 #. Post process (e.g. :meth:`smoothing <Table.smooth>`,
    :meth:`normalizing<Table.normalize>`),
    :meth:`plot the data <Table.write_asgl>`, or
    :meth:`write the table to a file <Table.write_hdf5>`.

A simple example, which simply collects the distribution of residue types in
a PDB file, is shown below:

.. literalinclude:: ../examples/residue_type.py
   :language: python

For more applied examples, see :ref:`study`.
