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

Utility functions
-----------------

.. autofunction:: uniform_bins

.. autofunction:: write_bondlib

.. autofunction:: write_anglelib

.. autofunction:: write_improperlib

.. autofunction:: write_splinelib

.. autofunction:: write_2dsplinelib

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
