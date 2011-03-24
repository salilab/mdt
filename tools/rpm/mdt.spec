Name:          mdt
Version:       SVN
Release:       1
License:       GPLv2
Summary:       Protein structure analysis
Group:         Applications/Engineering
Packager:      Ben Webb <ben@salilab.org>
Vendor:        Andrej Sali
URL:           http://www.salilab.org/mdt/
Source0:       %{name}-%{version}.tar.gz
BuildRoot:     %{_tmppath}/%{name}-%{version}-root
BuildRequires: modeller, glib2-devel >= 2.2.0
BuildRequires: python-devel, python-sphinx

%description
MDT prepares a raw frequency table, given information from MODELLER alignments
and/or PDB files. It can also process the raw frequency table in several ways
(e.g., normalization, smoothing, perform entropy calculations, and write out
the data in various formats, including for plotting by ASGL) and use as
restraints by MODELLER.

More precisely, MDT uses a sample of sequences, structures, and/or alignments
to construct a table N(a,b,c,...,d) for features a, b, c, ..., d. The sample
for generating the frequencies N is obtained depending on the type of features
a, b, c, ..., d. The sample can contain individual proteins, pairs of proteins,
pairs of residues in proteins, pairs of aligned residues, pairs of aligned
pairs of residues, chemical bonds, angles, dihedral angles, and pairs of
tuples of atoms. Some features work with triple alignments, too. All the
needed features a, b, c, ..., d are calculated automatically from the
sequences, alignments, and/or PDB files. The feature bins are defined
by the user when each feature is created.

%prep
%setup

%build
scons

%clean
[ "$RPM_BUILD_ROOT" != "/" ] && rm -rf ${RPM_BUILD_ROOT}

%check
scons test

%install
scons destdir=${RPM_BUILD_ROOT} docdir=/usr/share/doc/%{name}-%{version} \
      install docinstall

%files
%defattr(-,root,root)
%{_libdir}/libmdt.so
%{_libdir}/python2.*/site-packages/_mdt.so
%{_libdir}/python2.*/site-packages/mdt
%doc /usr/share/doc/%{name}-%{version}
/usr/share/mdt

%changelog
* Thu Mar 24 2011 Ben Webb <ben@salilab.org>   SVN-1
- First RPM build.
