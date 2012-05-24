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

libopenlat_la_SOURCES  = htk/lex.cpp htk/yacc.cpp htk-compiler.cpp iofilter.cpp
#libopenlat_la_LIBADD =  -lfst  -ldl -lz -lbz2 

htk/yacc.cpp htk/yacc.hpp htk/location.hh htk/position.hh htk/stack.hh: $(srcdir)/htk/htk.y
	$(YACC) $(YFLAGS) $(AM_YFLAGS) --output=htk/yacc.cpp --defines=htk/yacc.hpp $(srcdir)/htk/htk.y --report-file=htk/yacc.report 

htk/lex.o: htk/lex.cpp
	$(LIBTOOL) --mode=compile $(GCC)  $(CXXFLAGS) $(AM_CXXFLAGS) -c htk/lex.cpp

htk/lex.cpp: $(srcdir)/htk/htk.l
	$(LEX) $(LFLAGS) $(AM_LFLAGS) --outfile=htk/lex.cpp $(srcdir)/htk/htk.l 

