AM_CFLAGS = -std=c99 -I$(top_srcdir)/src/include
AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
#AM_LDFLAGS = -avoid-version -module -shared -export-dynamic -lltdl
#AM_LDFLAGS = -lltdl
AM_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB)

BUILT_SOURCES = htk/lex.cpp htk/yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh
CLEANFILES    = htk/lex.cpp htk/yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh htk/yacc.report

# yacc flags
AM_YFLAGS = --report=all 
# lex flags
AM_LFLAGS =

include_HEADERS = 

lib_LTLIBRARIES = libopenlat.la

libopenlat_la_SOURCES  = htk/lex.cpp htk/yacc.cpp htk-compiler.cpp iofilter.cpp sigsegv.cpp
#libopenlat_la_LIBADD =  -lfst  -ldl -lz -lbz2 

htk/yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh: $(srcdir)/htk/htk.y
	$(YACC) $(YFLAGS) $(AM_YFLAGS) --output=htk/yacc.cpp --defines=htk/yacc.hpp $(srcdir)/htk/htk.y --report-file=htk/yacc.report 

htk/lex.o: htk/lex.cpp
	$(LIBTOOL) --mode=compile $(GCC)  $(CXXFLAGS) $(AM_CXXFLAGS) -c htk/lex.cpp

htk/lex.cpp: $(srcdir)/htk/htk.l
	$(LEX) $(LFLAGS) $(AM_LFLAGS) --outfile=htk/lex.cpp $(srcdir)/htk/htk.l 


### Python extension ###
if HAVE_SWIG
BUILT_SOURCES = $(builddir)/openlat_swig_wrap.cpp
SWIG_SOURCES = openlat.i extra_typemaps.i

pkgpython_PYTHON = openlat.py
pkgpyexec_LTLIBRARIES   = _openlat.la
_openlat_la_SOURCES  = $(BUILT_SOURCES) $(SWIG_SOURCES)
_openlat_la_CXXFLAGS = $(COMPILANCE_CXXFLAGS) $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src/include 
_openlat_la_LDFLAGS  = -version-info 0:0:0 -export-dynamic -module $(PYTHON_LDFLAGS) $(LDFLAGS) -lfst

$(builddir)/openlat_swig_wrap.cpp : $(SWIG_SOURCES) 
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/src/include -o $@ $<
endif
