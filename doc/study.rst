.. highlight:: rest

.. currentmodule:: mdt

.. _study:

Sample studies with MDT
=======================

.. _introduction:

Introduction
------------

.. todo::

   Should have plots of raw data histograms superposed on the
   final restraints in all cases.

This section describes the use of MDT for updating many of the
MODELLER restraint libraries, including stereochemical, non-bonded, and
homology-derived restraints:

 #. Stereochemical restraints

    * chemical bonds: *p(Bond | BondType)*
    * chemical angles: *p(Angle | AngleType)*
    * improper dihedral angles as defined in the CHARMM residue
      topology file: *p(Dihedral | DihedralType)*
    * chemical angles: *p(Angle | AngleType)*
    * the ω dihedral angle of the peptide bond:
      *p(ω | ResidueType+1)* where ResidueType+1 refers to the residue type
      following the residue with the ω dihedral angle
    * the Φ and Ψ dihedral angles:
      *p(Φ | ResidueType)*, *p(Ψ | ResidueType)*
    * the sidechain dihedral angles:
      *p(χ1 | ResidueType)*, *p(χ2 | ResidueType)*, *p(χ3 | ResidueType)*,
      *p(χ4 | ResidueType)*
    * the mainchain conformation:
      *p(Φ, Ψ | ResidueType)*

 #. Non-bonded restraints

    * the mainchain hydrogen bonding restraints:
      *p(h | d, a)*
    * the non-bonded pair of atom triplets:
      *p(d, α1, α2, θ1, θ2, θ3 | t1, t2)*

 #. Homology-derived restraints

    * distance: *p(d | d')*

The following sections will outline the process of starting with the Protein
Data Bank (PDB) and ending up with the MODELLER restraint
library files. We will describe the rationale for the process, input data sets,
programs and scripts used, and even the analysis of the results. All of the
input files should be found in the MDT distribution, in the 
:file:`constr2005` directory.

The overall approach is to construct appropriately accurate, smooth, and
transformed surfaces based on the statistics in PDB for use as spatial
restraints for model building. The restraints from the first iteration will be
used to construct many models, which in turn will be used to re-derive the
equivalent restraints from the models. These model-derived restraints will
then be compared against the original PDB-derived restraints to find
problems and get indications as to how to change the restraints so that
models are statistically as similar to PDB structures as possible.

.. _stereo:

Stereochemical restraints
-------------------------

.. _sample:

The sample
~~~~~~~~~~

The starting point for deriving the restraints in this section consists of
9,365 chains that are representative of the 65,629 chains in the October 2005
version of PDB. The representative set was obtained by clustering all PDB
chains with MODELLER, such that the representative chains are from 30 to
3000 residues in length and are sharing less than 60% sequence identity to
each other (or are more than 30 residues different in length). This
is the corresponding MODELLER script:

.. literalinclude:: ../constr2005/cluster-PDB/make-pdb60.py
   :language: python

The actual chains for restraint derivation are in fact a subset of the 9,365
representative chains, corresponding to the 4,532 crystallographic structures
determined at least at 2 Å resolution (the representative structure for
each group is the highest resolution x-ray structure in the group). This
decision was made by looking at the distribution of the
χ1 dihedral angles as a function of resolution
(see :ref:`chi1`) and the distribution of resolutions themselves
for all 9,365 representative chains, using this MDT script:

.. literalinclude:: ../constr2005/impact-of-resolution/make-mdt.py
   :language: python

This script creates a :class:`Library` object and then adds an X-ray
resolution feature. Values of this feature are placed into 20 bins of width 0.2,
starting at 0. It then creates a
:class:`Table` object, which is the MDT table itself.
This starts off as an empty 1D table of the X-ray
resolution feature. It then uses a MODELLER alignment object to read the
sequences from :file:`pdb_60.pir` one by one, and for each one it
updates the X-ray resolution feature in the MDT table by calling
:meth:`Table.add_alignment`. Finally,
the table is written out to the file :file:`mdt2.mdt` using
:meth:`Table.write`.

The resulting MDT table :file:`mdt2.mdt` was then
plotted with the script:

.. literalinclude:: ../constr2005/impact-of-resolution/asgl.py
   :language: python

where the :meth:`Table.write_asgl`
method writes out an ASGL script and the MDT data in a form suitable for
plotting (which we then execute with ASGL using Python's
:meth:`os.system` method). This gives an
:download:`impact of resolution plot <../constr2005/impact-of-resolution/asgl2-a.pdf>`.


.. _chembonds:

Chemical bonds
~~~~~~~~~~~~~~

The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/bonds/make-mdt.py
   :language: python

In this case, we read the file :file:`bndgrp.lib` which defines all
chemical bonds, using the :meth:`BondClasses.read` method. The
MDT we then construct is a 3D table of X-ray resolution, bond type, and bond
length. The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/bonds/asgl.py
   :language: python

giving :download:`a set of bond plots <../constr2005/bonds/asgl1-a.pdf>`.
Notice that here we use the :meth:`Table.reshape` method,
which can reshape a table by reordering the features, and/or reducing the bin
ranges (offset or shape) of these features. In this case we don't change the
feature order, or the offset (leaving it at the default 0,0,0) but we do
change the shape. The first feature is restricted to only one bin - because
our X-ray resolution feature contains only two bins (for "less than
2 Å" and the undefined bin, which catches everything 2 Å or
greater) this keeps only the good structures
for our plot. The other two features have their bin ranges reduced by 1
(a negative value for shape means "reduce the size by this amount"),
which effectively removes the final ("undefined") bin.

Inspection of the plots shows that all distributions are mono-modal, but most
are distinctly non-Gaussian. However, at this point, Gaussian restraints are
still favored because the ranges are very narrow and because the non-Gaussian
shape of the histograms may result from the application of all the other
restraints (this supposition will be tested by deriving the corresponding
distributions from the models, not PDB structures). Also, although many
distributions are quite symmetrical, not all of them are. Therefore, there is
the question of how best to fit a restraint to the data. There are at least
three possibilities, in principle: (i) calculating the average and standard
deviation from all (subset) of the data, (ii) least-squares fitting of the
Gaussian model to the data, and (iii) using cubic splines of the data. The
first option was adopted here: the mean and standard deviation will be the
parameters of the analytically defined bond restraint for MODELLER.

The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/bonds/modlib.py
   :language: python

Here we use the :meth:`Table.integrate`
method, which removes a feature from the table by integrating the remaining
features over all of that feature's bins, and the
:func:`write_bondlib` function
to write out a MODELLER script which builds restraints using the MDT-derived
distributions.

.. _chemangles:

Chemical angles
~~~~~~~~~~~~~~~

As for the bonds above, the MDT table is constructed with the following
MDT Python script:

.. literalinclude:: ../constr2005/angles/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/angles/asgl.py
   :language: python

giving :download:`a set of angle plots <../constr2005/angles/asgl1-a.pdf>`.

The situation is similar to that for the chemical bonds, except that there are
also four cases of bi-modal (as opposed to mono-modal) distributions:
Asp:CB:CG:OD2, Asp:OD2:CG,OD1, Pro:CB:CG:CD, and Pro:CD:N:CA angles. The Asp
bi-modal distribution may result from crystallographers mislabeling carboxyl
oxygens for the protonated state of the sidechain (which is interesting
because the corresponding angles might be used as a means to assign the
protonation state). The mean values for these angles were edited by hand.
Otherwise exactly the same considerations as for bonds apply here.

The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/angles/modlib.py
   :language: python

.. _impropers:

Improper dihedral angles
~~~~~~~~~~~~~~~~~~~~~~~~

Exactly the same situation applies as for the chemical bonds.
The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/impropers/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/impropers/asgl.py
   :language: python

giving
:download:`a set of improper plots <../constr2005/impropers/asgl1-a.pdf>`.

The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/impropers/modlib.py
   :language: python

.. _chi1:

Sidechain dihedral angle χ1
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first question asked was "What is the impact of resolution on the
distributions of residue χ1?". The answer was
obtained by constructing and inspecting
*p(χ1 | R, resolution)* with:

.. literalinclude:: ../constr2005/chi1/impact-of-resolution/make-mdt.py
   :language: python

and

.. literalinclude:: ../constr2005/chi1/impact-of-resolution/asgl.py
   :language: python

giving
:download:`this output <../constr2005/chi1/impact-of-resolution/asgl2-a.pdf>`
which clearly shows that X-ray structures at resolution of at least 2.0 Å
are just fine. X-ray structures above that resolution and NMR structures
(whose resolution is set artificially to 0.45 Å for the purposes of MDT
tabulation) do not appear to be suitable for deriving restraints for modeling,
as the peaks are significantly wider and there is a significant population
at ~120°. Also, the peaks appear Gaussian. Thus, a weighted sum of three
Gaussians (except for Pro, which has two) was judged to be an appropriate
model for these data. Again, the following script was used to construct the
MDT table:

.. literalinclude:: ../constr2005/chi1/make-mdt.py
   :language: python

and the contents then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi1/asgl.py
   :language: python

giving :download:`a set of χ1 plots <../constr2005/chi1/asgl1-a.pdf>`.

The weights, means, and standard deviations of the Gaussians were obtained
by least-squares fitting with ASGL (with the script below) and are manually
added to the MODELLER MDT library.

.. literalinclude:: ../constr2005/chi1/fit.top

.. _chi2:

Sidechain dihedral angle χ2
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The situation is very similar to that for χ1,
except that the shapes of histograms are not Gaussian in most cases.
Therefore, 1D cubic splines are used to represent the restraints.

The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/chi2/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi2/asgl.py
   :language: python

giving :download:`a set of χ2 plots <../constr2005/chi2/asgl1-a.pdf>`.

The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/chi2/modlib.py
   :language: python

This script also uses :meth:`Table.smooth` to smooth the raw
frequencies and :meth:`Table.normalize` to convert the
distribution into a PDF. It is then converted into a statistical potential
by taking the negative log of the values (using the
:meth:`Table.log_transform`,
:meth:`Table.linear_transform`, and
:meth:`Table.offset_min` methods).
The smoothing parameter *weight* of 10 was selected by trial and error,
inspecting the resulting restraints in :file:`modlib-a.ps`, also produced
by the script above.

.. _chi3:

Sidechain dihedral angle χ3
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exactly the same considerations apply as to χ2.
The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/chi3/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi3/asgl.py
   :language: python

giving :download:`a set of χ3 plots <../constr2005/chi3/asgl1-a.pdf>`.
The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/chi3/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _chi4:

Sidechain dihedral angle χ4
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exactly the same considerations apply as to χ2 and χ3. The MDT table is
constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/chi4/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi4/asgl.py
   :language: python

giving :download:`a set of χ4 plots <../constr2005/chi4/asgl1-a.pdf>`.
The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/chi4/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _phi:

Mainchain dihedral angle Φ
~~~~~~~~~~~~~~~~~~~~~~~~~~

Exactly the same considerations apply as to χ2, χ3, and χ4. The MDT
table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/phi/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/phi/asgl.py
   :language: python

giving :download:`a set of Φ plots <../constr2005/phi/asgl1-a.pdf>`.
The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/phi/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _psi:

Mainchain dihedral angle Ψ
~~~~~~~~~~~~~~~~~~~~~~~~~~

Exactly the same considerations apply as to χ2, χ3, χ4, and Φ.
The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/psi/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/psi/asgl.py
   :language: python

giving :download:`a set of Ψ plots <../constr2005/psi/asgl1-a.pdf>`.
The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/psi/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _omega:

Mainchain dihedral angle ω
~~~~~~~~~~~~~~~~~~~~~~~~~~

This dihedral angle is a little different from all others explored thus far
because it depends more strongly on the type of the subsequent residue than
the type of the residue whose dihedral angle is studied; that is, the ω
of the residue preceding Pro, not the Pro ω, is impacted by the Pro
residue. These dependencies are explored with MDT tables in directory
:file:`constr2005/omega/run1/`. The bottom line is that we need
to set *delta* to 1 when creating our
:class:`residue type <features.ResidueType>` feature (rather
than the default value 0), which will make it refer to the type of the
residue after the residue with the dihedral angle ω.

The next step is to obtain the *p(ω | R+1)*
distributions with finer sampling of 0.5°:

.. literalinclude:: ../constr2005/omega/make-mdt.py
   :language: python

The :download:`distribution in raw form <../constr2005/omega/asgl1-a.pdf>`
is then plotted with:

.. literalinclude:: ../constr2005/omega/asgl.py
   :language: python

and :download:`in logarithmic form <../constr2005/omega/asgl2-a.pdf>`
with:

.. literalinclude:: ../constr2005/omega/asgl-log.py
   :language: python

Clearly, the peaks are sharp and will best be modeled by Gaussian distributions.

Similarly to χ1, two Gaussian distributions are
fit to the histograms with the following ASGL script:

.. literalinclude:: ../constr2005/omega/fit.top

The means and standard deviations for each residue type are extracted from
:file:`fit.log` by the ASGL :file:`get_prms.F`
program, but they are only used to guess the means of 179.8° and 0°
and standard deviations of 1.5° and 1.5° for the two peaks,
respectively. The distribution of ω dihedral angles in the models
calculated with these ω restraints will be checked carefully and the
restraint parameters will be adapted as needed.  

The weights of the peaks are not determined reliably by least-squares fitting
in this case because the second weight is very close to 0 (in principle, they
can even be less than zero). Therefore, they are determined separately by
establishing *p(cω | R+1)* where *cω* is the class of the
ω dihedral angle (1 or 2, *trans* or *cis*).

The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/omega/weights/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/omega/weights/asgl.py
   :language: python

giving
:download:`an omega weights plot <../constr2005/omega/weights/asgl1-a.pdf>`.

The library :file:`omega.py` is edited manually to replace the
means and standard deviations with `179.8  0.0   2.3  2.3`.

.. _phipsi:

Mainchain dihedral angles Φ and Ψ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initial runs in :file:`run1` explored Ramachandran maps
extracted from different representative sets of structures (e.g., clustered by
40% sequence identity) and stratification by the crystallographic residue
B\ :sub:`iso` as well as resolution and residue type. We ended up with
the sample and stratification described above.

The 2D histograms *p(Φ, Ψ | R)* are derived with:

.. literalinclude:: ../constr2005/phipsi/make-mdt.py
   :language: python

They are plotted with

.. literalinclude:: ../constr2005/phipsi/asgl.py
   :language: python

giving :download:`a set of Φ/Ψ plots <../constr2005/phipsi/asgl1-a.pdf>`.

The distributions are clearly not 2D Gaussian functions and need to be
approximated by 2D cubic splines. Exploring and visualizing various smoothing
options results in the following file to produce the final MODELLER MDT
library:

.. literalinclude:: ../constr2005/phipsi/modlib.py
   :language: python

The raw, smooth, and transformed surfaces are visualized and compared best
with Mathematica.

.. _nonbonded:

Non-bonded restraints
---------------------

A general pairwise distance- and atom-type dependent statistical potential
for all atom type pairs has been derived by Min-yi Shen (DOPE).
MDT could, however, be used to derive specialized pairwise non-bonded
restraints.
