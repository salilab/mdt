%module _mdt
%include "typemaps.i"

%{
#include "../src/mdt.h"
%}

%apply int *OUTPUT { int *ierr };

%include "../src/mdt.h"
