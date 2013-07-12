%module(directors="1") openlat 

%{
#define SWIG_FILE_WITH_INIT
#include <openlat/normalize.h>
using namespace openlat;
using namespace fst;

%}



%include <typemaps.i>
%include <stl.i>
%include <std_pair.i>
%include <std_string.i>
%include <std_vector.i>
%include "extra_typemaps.i"

using namespace std;


// the directors must be parsed before they are used by others
%include <openlat/normalize.h>




/*****************************************************************
 * Macro for instantiating templates for all the different classes
 * and functions.  Note that you must arrange for
 *
 * to be explicitly defined beforehand.
 *****************************************************************/

%define INSTANTIATE_FST(ARCTYPE)

void Normalize(fst::MutableFst<ARCTYPE> *fst, float posterior_scale = 1.0) ;
float Entropy(const fst::Fst<ARCTYPE> &fst);

%enddef


INSTANTIATE_FST(StdArc);
INSTANTIATE_FST(LogArc);
