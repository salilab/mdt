.. highlight:: rest

The :mod:`mdt` Python module
================================

.. automodule:: mdt

.. autoclass:: mdt.Library
   :members: atom_classes, bond_classes, angle_classes, dihedral_classes,
             tuple_classes, hbond_classes

.. autoclass:: mdt.Table

.. autoclass:: mdt.TupleClasses

.. autoclass:: mdt.BondClasses

.. autoclass:: mdt.HydrogenBondClasses

.. autoclass:: mdt.TableSection

.. autofunction:: uniform_bins

.. autofunction:: write_anglelib

.. autofunction:: write_improperlib

.. autofunction:: write_splinelib

.. autofunction:: write_2dsplinelib

.. data:: Float

.. exception:: mdt.MDTError

   A generic MDT error.

.. exception:: mdt.FileFormatError

   File format error.
