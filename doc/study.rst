.. highlight:: rest

.. currentmodule:: mdt

.. _study:

Sample studies with MDT
=======================

.. _introduction:

Introduction
------------

.. note::

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
=========================

.. _sample:

The sample
----------

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
`impact of resolution plot <pdf/impact-of-resolution.pdf>`_.


.. _chembonds:

Chemical bonds
--------------

The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/bonds/make-mdt.py
   :language: python

In this case, we read the file :file:`bndgrp.lib` which defines all
chemical bonds, using the :meth:`BondClasses.read` method. The
MDT we then construct is a 3D table of X-ray resolution, bond type, and bond
length. The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/bonds/asgl.py
   :language: python

giving `a set of bond plots <pdf/bonds.pdf>`_. Notice that
here we use the :meth:`Table.reshape` method,
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
---------------

As for the bonds above, the MDT table is constructed with the following
MDT Python script:

.. literalinclude:: ../constr2005/angles/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/angles/asgl.py
   :language: python

giving `a set of angle plots <pdf/angles.pdf>`_.

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
------------------------

Exactly the same situation applies as for the chemical bonds.
The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/impropers/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/impropers/asgl.py
   :language: python

giving `a set of improper plots <pdf/impropers.pdf>`_.

The final MODELLER MDT library is produced with:

.. literalinclude:: ../constr2005/impropers/modlib.py
   :language: python

.. _chi1:

Sidechain dihedral angle χ1
---------------------------

The first question asked was "What is the impact of resolution on the
distributions of residue χ1?". The answer was
obtained by constructing and inspecting
*p(χ1 | R, resolution)* with:

.. literalinclude:: ../constr2005/chi1/impact-of-resolution/make-mdt.py
   :language: python

and

.. literalinclude:: ../constr2005/chi1/impact-of-resolution/asgl.py
   :language: python

giving `this output <pdf/chi1-impact.pdf>`_
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

giving `a set of χ1 plots <pdf/chi1.pdf>`_.

The weights, means, and standard deviations of the Gaussians were obtained
by least-squares fitting with ASGL (with the script below) and are manually
added to the MODELLER MDT library.

.. literalinclude:: ../constr2005/chi1/fit.top

.. _chi2:

Sidechain dihedral angle χ2
---------------------------

The situation is very similar to that for χ1,
except that the shapes of histograms are not Gaussian in most cases.
Therefore, 1D cubic splines are used to represent the restraints.

The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/chi2/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi2/asgl.py
   :language: python

giving `a set of χ2 plots <pdf/chi2.pdf>`_.

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
---------------------------

Exactly the same considerations apply as to χ2.
The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/chi3/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi3/asgl.py
   :language: python

giving `a set of χ3 plots <pdf/chi3.pdf>`_. The final MODELLER MDT library is
produced with:

.. literalinclude:: ../constr2005/chi3/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _chi4:

Sidechain dihedral angle χ4
---------------------------

Exactly the same considerations apply as to χ2 and χ3. The MDT table is
constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/chi4/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/chi4/asgl.py
   :language: python

giving `a set of χ4 plots <pdf/chi4.pdf>`_. The final MODELLER MDT library is
produced with:

.. literalinclude:: ../constr2005/chi4/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _phi:

Mainchain dihedral angle Φ
--------------------------

Exactly the same considerations apply as to χ2, χ3, and χ4. The MDT
table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/phi/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/phi/asgl.py
   :language: python

giving `a set of Φ plots <pdf/phi.pdf>`_. The final MODELLER MDT library
is produced with:

.. literalinclude:: ../constr2005/phi/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _psi:

Mainchain dihedral angle Ψ
--------------------------

Exactly the same considerations apply as to χ2, χ3, χ4, and Φ.
The MDT table is constructed with the following MDT Python script:

.. literalinclude:: ../constr2005/psi/make-mdt.py
   :language: python

The contents of the MDT table are then plotted with ASGL as follows:

.. literalinclude:: ../constr2005/psi/asgl.py
   :language: python

giving `a set of Ψ plots <pdf/psi.pdf>`_. The final MODELLER MDT library
is produced with:

.. literalinclude:: ../constr2005/psi/modlib.py
   :language: python

The resulting restraints are plotted in :file:`modlib-a.ps`,
also produced by the script above.

.. _omega:

Mainchain dihedral angle ω
--------------------------

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

The `distribution in raw form <pdf/omega.pdf>`_ is then plotted with:

.. literalinclude:: ../constr2005/omega/asgl.py
   :language: python

and `in logarithmic form <pdf/omega-log.pdf>`_ with:

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

giving `an omega weights plot <pdf/omega-weights.pdf>`_.

The library :file:`omega.py` is edited manually to replace the
means and standard deviations with `179.8  0.0   2.3  2.3`.

.. _phipsi:

Mainchain dihedral angles Φ and Ψ
---------------------------------

The initial runs in :file:`run1` explored Ramachandran maps
extracted from different representative sets of structures (e.g., clustered by
40% sequence identity) and stratification by the crystallographic residue
Biso as well as resolution and residue type. We ended up with the sample
and stratification described above.

The 2D histograms *p(Φ, Ψ | R)* are derived with:

.. literalinclude:: ../constr2005/phipsi/make-mdt.py
   :language: python

They are plotted with

.. literalinclude:: ../constr2005/phipsi/asgl.py
   :language: python

giving `a set of Φ/Ψ plots <pdf/phipsi.pdf>`_.

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
=====================

A general pairwise distance- and atom-type dependent statistical potential
for all atom type pairs has been derived by Min-yi Shen (DOPE). Here, we focus
on specialized pairwise non-bonded restraints.

.. _hbond:

Mainchain hydrogen bonding restraints
-------------------------------------

The idea is to describe them as restraints on the donor and acceptor pairs of
atom triplets. The donor triplet could be N\ :sub:`i` -
CA\ :sub:`i` - C\ :sub:`i-1` and the acceptor
triplet could be O\ :sub:`j` - C\ :sub:`j` -
CA\ :sub:`j`, where i and j are residue indices.

As always, here are the aspects that need to be explored and defined:

Which dependent features to use?
   The dependent features are clearly the distance, two angles, and three
   dihedral angles of the two triplets of atoms (the donor and acceptor
   triplets), though I hope that some of them can be omitted without too much
   of a problem.

   We could have the total potential as a sum of terms for each one of the
   dependent features, but then the correlations between them would be lost.
   We need to find out which features are most correlated and join those in the
   higher order restraints.

   Physically, how does the definition of the H atom position on the donor
   'eliminate' the need to consider the two atoms connected to the the donor N
   (and thus reduce the number of dependent features)?

Which independent features to use?
   The independent features can be divided into those that the dependent
   features 'really' depend on and those that are there for quality control
   (e.g., resolution). The independent features include the triplet types
   (donor and acceptor, irrespective of the residue type), sequence separation,
   and X-ray structure resolution. It seems best to fix the triplet types to DON
   and ACC respectively and let the sequence separation span the negative and
   positive range. This way, the triplet types could even be omitted from
   the MDT table (they don't change), just like the resolution (though I
   don't like the former omission).

What is the range and binning of these features?
   The first (4.5 Å) and second atom shell (8 Å) are important
   numbers for considering the range of the dependent distance. In addition,
   the standard range of H-bonds is 3.5 Å.

   It is conceivable that looking only at the raw distributions for deciding
   about the dependent and independent features, their bins, and ranges would
   be misleading. For example, normalization of the raw frequency with that
   expected by chance might eliminate a large number of differences caused by
   such features as sequence separation. Thus, it may be possible to use
   coarser binning in some independent features.

   The previous point indicates the need to develop smoothing and normalization
   early on, so that 'final' restraints and not raw frequencies are used in
   judging the selection of features, binning, and range. One should be helped
   by conditional entropies and Mathematica in this endeavor.

How to smooth the raw frequencies?
   By adding a uniform distribution with an appropriate weight.

How to normalize the raw frequencies?
   The problem is that both the 'analytical' and 'empirical' routes are very
   difficult: (i) duplicating the 4 π r\ :sup:`2` argument here would require
   considering volume elements spanned by distance, dihedral angles, and
   angles, which is difficult; (ii) it is difficult to imagine what pairs
   of atom triplets in real structures would provide a good reference.
   Min-yi probably came to the rescue with an idea to simulate
   pairs of triplets inside a sphere of say 23 Å radius. Pairs of triplets
   are placed randomly inside the sphere, no atom-overlap checks are performed,
   and then the distribution of the relative orientations is collected. But
   there is still a problem with this idea: Because the reference does not
   depend on sequence separation, the reference will not 'normalize' out the
   impact of sequence separation in the raw frequencies, which does not
   'feel' right.

   Here is the beginning of a larger idea, based on playing games with a system
   we define, so we know what it is. It is also based on an idea of progression,
   evolving the system from a simple version where everything is clear to a more
   complex and more realistic system in a series of steps that are hopefully
   managable. And it is already clear up front that the idea will be fighting
   both the thermodynamic assumption for deriving the statistical potentials
   and the multi-body problem (because the idea is principled and because these
   are the two main issues in extracting restraints from a sample of
   structures).

   So the starting toy system is a polypeptide chain that 'feels' energy terms
   for only chemical bonds, angles, and dihedral angles, each one of which
   depends on the atom types and the residue type. There are no other
   interactions, not even non-bonded interactions. The chain looks like a
   random walk (but of course it is not). By definition, the 'PDB' (i.e., the
   sample) contains native structures at the global energy minimum, each one of
   which is entirely self-consistent with each other (i.e., there is no
   frustration among the restraints). Clearly, the sample will show that the
   distributions of the bonds, angles, and dihedral angles are delta functions.
   Therefore, we can in a straighforward way determine the means of all
   restraints, but there is no 'entropy' in the resulting restraints and/or it
   is not knowable from the sample. Nevertheless, the corresponding pdf's (i.e.,
   delta functions) would allow us to exactly predict the native
   structure of any new sequence.

   But wait a second, we just may have tacitly skipped the 'normalization' step
   because the final answer was so obvious. Should we in fact formally normalize
   the delta functions by a distribution of bonds, angles, and dihedral angles
   for a random collection of points (of course, we would get the same delta
   functions back)? Why a random collection? Because it seems that the reference
   distribution should be based on all energy terms but those we are trying to
   extract from the sample. To make these steps clearer, let's consider them
   only at the 2nd level of buildup, next. But we should come back here and
   define exactly the properties of the random collection (e.g., What is random?
   Uniformly distributed in Cartesian coordinates? What is volume shape and
   size? These will impact on the distance distribution, though not on the angle
   and dihedral angle distributions if the volume is large enough.).

   Is this progression very similar to Min-yi's DOPE manuscript, except that
   partitioning is a little different (there, it is from single-body to
   two-body)?

   Now comes the step to the 2nd level of the buildup. Let's say that the real
   chains also feel the non-bonded Lennard-Jones terms between all atoms
   separated by more than 3 chemical bonds, in addition to the bond, angle, and
   dihedral angle terms. The systems is now frustrated and determining the
   Lennard-Jones terms is not easy (if we had hard spheres, then we could just
   look for the closest distance in a large sample and be done). Why do I feel
   that -k\ :sub:`B`\ T times the natural logarithm of the ratio
   between the PDB sample distribution of non-bonded distances and the
   non-bonded distribution from the sample in the first step is a good
   approximation to the Lennard-Jones energy (not even PMF, but real potential
   energy)?

   In the end, the order of the steps seems important for the final result,
   while physically it should not be. Should we iterate or explore multiple
   paths for a self-consistent solution? It seems estimating the more
   determining terms (smaller entropy on their own) first makes sense.

   Do it for homology-derived restraints, which are not physical, but the
   statistical framework above should still apply. They are strong restraints,
   so they should probably be considered early in the progression, right after
   the stereochemistry.

   #. For an ideal system in which all variables are independent
      (uncoupled), each variable will be found in its minimum/equilibrium
      position. No information about force constants will be revealed.

   #. Introducing an distance restraint into a protein chain will
      inevitably include a variable coupled to the existing internal
      coordinates into the minimization, i.e. the restraint "bond length" is
      not an independent variable.

   #. The addition of the restraints will expose the force constants for
      each existing restraint.

   Then some comments on your writeup, which is very cool! I like the specific
   example of the simplest possible example of a "frustrating" restraint (i.e.,
   the restraint on l), with equations. Did you try Mathematica to find the
   solution of the system of three equations? I am sure they can be solved
   numerically at the very least. There are probably no insights in the specific
   solution, which must depend, qualitatively even,  on the values of the means
   and force constants, but who knows. Also, one could write other systems of
   three equations by picking say distance l, not theta, as the variable to
   minimize. I wonder if that would give us any advantage?

   * I certainly agree completely with your 3 points below. it is interesting
     that this "gedanken" experiment addresses both the multiple body
     problem as well as the concern that PDB is not a boltzmann ensemble.
     This leads me to believe that we are doing something fundamentally correct.

   * It seems you can improve the description of the reference state and
     normalization in the dope paper, based on the discussion we have here.
     Again, I think it would be best to keep the motivation/rationale/execution
     of normalization as statistical as possible, as opposed to physical.
     Or at least do so at the beginning, and then make the physical
     connection in the end, if you must, just like you did overall.
     On the other hand, our discussion here is physical, not statistical,
     so i am not sure about the comment i just made.

   * Still for the DOPE paper, specifically: why did we select the reference
     state for DOPE the way we did (which interactions are on/off ...).

   * I wonder if going from 3 to 4 points would again fundamentally change the
     situation in 3D (as it did when we went from 2 to 3 points), since we would
     for the first time introduce chirality. it probably does not matter in ways
     relevant to us here.

   * So now to assumptions/approximations/new questions (for you ;-) ):
     I suppose we won't be able to solve the problem exactly (actually,
     it would be very good to exactly define the problem we are addressing:
     How do we extract the most accurate means and force constants from
     a sample of native states in our "toy" universe?). So we need to
     start thinking about suitable assumptions/approximations. What
     follows uses language so unprecise it is irresponsible to use it.
     But what can i do! ;-)

   * It may be possible to come up with reasonable approximations if we assume
     that the force constants for the added "frustrating" restraints
     (e.g., k\ :sub:`l`) are very weak compared to the other
     restraints (e.g., bonds, angles, dihedrals).

   * So maybe in our gedanken experiment we can add increasingly weaker
     restraints and extract them from the native states and reference
     states corresponding to global minima from the previous step (what
     we had on the whiteboard in my office). in addition, maybe a "mean
     field" solution is the best we can do; ie, after all unsolvable
     minimum-defining equations end up pushing the means and force
     constants in all directions, we end up with some kind of a gaussian
     distribution for them and maybe we CAN estimate the mean
     and standard deviation of that gaussian, although we cannot get the
     individual values in it.

   * So I'd like to ask you here again how come people use E = -kT
     ln p(native) / p(reference) ; what exactly is the native, reference, E for
     which this equation allows one to calculate E. what is the origin of this
     equation, approximations, ...?

   * It still seems good to derive the expression equivalent to E = -kT
     ln p(native) / p(reference) explicitly for your simple 3 body system. Let's
     push that one to its complete/clear/explicit solution to improve our
     understanding in general. But you will need to make it more complicated by
     introducing a few different types of points, you need to imagine a PDB
     native triangle structures for lots of triangles consisting of the minimal
     number of types of vertices (just enough to make it useful here). Then,
     how do we get all the means and force constants?

Overall considerations
   There is a picture that a short-range, residue type-independent H-bonding
   potential will be useful for restraining the specific local geometry and
   not for selecting the fold directly. A longer range orientation-dependent
   two-body term could be helpful in discriminating between different folds.

   Presumably, the final restraints will not have a relatively small number of
   simple (Gaussian) peaks, thus the new XOR restraint is probably not
   indicated here.

   Develop some knowledge about the problem by doing lots of 1D histograms
   for the individual features in your sample, both dependent and independent
   ones.

.. _odnonbond:

Orientation-dependent pairwise non-bonded restraints
====================================================

.. _homology:

Homology-derived restraints
===========================
