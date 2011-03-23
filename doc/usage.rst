.. highlight:: rest

.. currentmodule:: mdt

Usage
=====

MDT is simply a Python extension module, and as such can be used in
combination with other Python modules, such as MODELLER or the Python
standard library.

Running pre-built binaries
--------------------------

In the Sali lab, MDT is built as part of the nightly build system, and is
installed with MODELLER in `/diva1/home/modeller`. Thus, you can run an MDT
script :file:`foo.py` just like any regular Python script with a command
similar to::

   /diva1/home/modeller/modpy.sh python foo.py

.. _compilation:

Compilation from source code
----------------------------

You can get the current MDT code by running the following::

   svn co https://svn.salilab.org/impmod/trunk impmod

If you already have a copy of MDT, you can update it to the current code
by running::

   svn update

To compile, create a file called :file:`config.py` in the
:file:`impmod/mdt` directory, and in it set the
*modeller* Python variable to the directory where you have
MODELLER installed. Then run `scons` in the same directory
(and optionally `scons test`) to build MDT. This will produce a script
:file:`bin/mdtpy.sh` which can be used to run an MDT Python script
:file:`foo.py`::

   bin/mdtpy.sh python foo.py

If you want to install MDT, run `scons install`. You can additionally specify
a `prefix` option to install in a different directory. For example,
`scons prefix=/foo install` will install MDT in the :file:`/foo` directory.

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
