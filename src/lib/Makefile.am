AM_CFLAGS = -std=c99 -I$(top_srcdir)/src/include
AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
#AM_LDFLAGS = -avoid-version -module -shared -export-dynamic -lltdl
#AM_LDFLAGS = -lltdl
AM_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB)

BUILT_SOURCES  = htk/htk-lex.cpp htk/htk-yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh
CLEANFILES     = htk/htk-lex.cpp htk/htk-yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh htk/yacc.report
#BUILT_SOURCES += thot/thot-lex.cpp thot/thot-yacc.cpp thot/yacc.hpp thot/location.hh thot/position.hh thot/stack.hh
#CLEANFILES    += thot/thot-lex.cpp thot/thot-yacc.cpp thot/yacc.hpp thot/location.hh thot/position.hh thot/stack.hh thot/yacc.report

# yacc flags
AM_YFLAGS = --report=all 
# lex flags
AM_LFLAGS =

include_HEADERS = 

lib_LTLIBRARIES = libopenlat.la

libopenlat_la_SOURCES  = htk/htk-lex.cpp htk/htk-yacc.cpp htk-compiler.cpp iofilter.cpp sigsegv.cpp
if HAVE_THOT
libopenlat_la_SOURCES  += thot-compiler.cpp
endif
#libopenlat_la_LIBADD =  -lfst  -ldl -lz -lbz2 

# htk parser and compiler
htk/htk-yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh: $(srcdir)/htk/htk.y
	$(YACC) $(YFLAGS) $(AM_YFLAGS) --output=htk/htk-yacc.cpp --defines=htk/yacc.hpp $(srcdir)/htk/htk.y --report-file=htk/yacc.report # -DBISON_VERSION=$(BISON_VERSION)

htk/htk-lex.o: htk/htk-lex.cpp
	$(LIBTOOL) --mode=compile $(GCC)  $(CXXFLAGS) $(AM_CXXFLAGS) -c htk/htk-lex.cpp

htk/htk-lex.cpp: $(srcdir)/htk/htk.l
	$(LEX) $(LFLAGS) $(AM_LFLAGS) --outfile=htk/htk-lex.cpp $(srcdir)/htk/htk.l 


## thot parser and compiler
#thot/thot-yacc.cpp thot/yacc.hpp thot/location.hh thot/position.hh thot/stack.hh: $(srcdir)/thot/thot.y
#	$(YACC) $(YFLAGS) $(AM_YFLAGS) --output=thot/thot-yacc.cpp --defines=thot/yacc.hpp $(srcdir)/thot/thot.y --report-file=thot/yacc.report 
#
#thot/thot-lex.o: thot/thot-lex.cpp
#	$(LIBTOOL) --mode=compile $(GCC)  $(CXXFLAGS) $(AM_CXXFLAGS) -c thot/thot-lex.cpp
#
#thot/thot-lex.cpp: $(srcdir)/thot/thot.l
#	$(LEX) $(LFLAGS) $(AM_LFLAGS) --outfile=thot/thot-lex.cpp $(srcdir)/thot/thot.l 


### Python extension ###
if HAVE_SWIG
BUILT_SOURCES += openlat_swig_wrap.cpp
SWIG_SOURCES = openlat.i extra_typemaps.i

pkgpython_PYTHON = openlat.py
pkgpyexec_LTLIBRARIES   = _openlat.la
_openlat_la_SOURCES  = openlat_swig_wrap.cpp $(SWIG_SOURCES)
_openlat_la_CXXFLAGS = $(COMPILANCE_CXXFLAGS) $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src/include 
_openlat_la_LDFLAGS  = -version-info 0:0:0 -export-dynamic -module $(PYTHON_LDFLAGS) $(LDFLAGS) -lfst

openlat_swig_wrap.cpp : $(SWIG_SOURCES) 
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/src/include -o $@ $<
endif
