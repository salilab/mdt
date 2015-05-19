.. highlight:: rest

The :mod:`mdt` Python module
================================

.. automodule:: mdt

Setup of the MDT system
-----------------------

.. autoclass:: Library
   :members:

.. autoclass:: TupleClasses
   :members:

.. autoclass:: BondClasses
   :members:

.. autoclass:: HydrogenBondClasses
   :members:

Creation and manipulation of data tables
----------------------------------------

.. autoclass:: Table
   :members:

.. autoclass:: TableSection
   :members:

.. autoclass:: Feature
   :members:

.. autoclass:: Bin
   :members:

.. autoclass:: Source
   :members:

Library information
-------------------

.. data:: version

   The full MDT version number, as a string, e.g. '5.0' or 'SVN'.

.. data:: version_info

   For release builds, the major and minor version numbers as a tuple of
   integers - e.g. (5, 0). For SVN builds, this is the same as 'version'.

Utility functions
-----------------

.. autofunction:: uniform_bins

.. autofunction:: write_bondlib

.. autofunction:: write_anglelib

.. autofunction:: write_improperlib

.. autofunction:: write_splinelib

.. autofunction:: write_2dsplinelib

.. autofunction:: write_statpot

Bin storage types
-----------------

.. data:: Float
          Double
          Int32
          UnsignedInt32
          Int16
          UnsignedInt16
          Int8
          UnsignedInt8

          See :ref:`binstorage`.

Exceptions
----------

.. exception:: MDTError

   A generic MDT error.

.. exception:: FileFormatError

   A file is of the wrong format.
