AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
#AM_LDFLAGS = -avoid-version -module -shared -export-dynamic -lltdl
#AM_LDFLAGS = -lltdl
AM_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB)

BUILT_SOURCES = lex.htk.cpp yacc.htk.cpp yacc.htk.hpp
CLEANFILES    = lex.htk.cpp yacc.htk.cpp yacc.htk.hpp yacc.htk.report

# yacc flags
AM_YFLAGS = --report=all 
# lex flags
AM_LFLAGS =

include_HEADERS = 

lib_LTLIBRARIES = libopenlat.la

libopenlat_la_SOURCES  = lex.htk.cpp yacc.htk.cpp htk-compiler.cpp iofilter.cpp
#libopenlat_la_LIBADD =  -lfst  -ldl -lz -lbz2 

yacc.htk.cpp yacc.htk.hpp: $(srcdir)/htk.y
	$(YACC) $(YFLAGS) $(AM_YFLAGS) --output=yacc.htk.cpp --defines=yacc.htk.hpp $(srcdir)/htk.y --report-file=yacc.htk.report 

lex.htk.o: lex.htk.cpp
	$(LIBTOOL) --mode=compile $(GCC)  $(CXXFLAGS) $(AM_CXXFLAGS) -c lex.htk.cpp

lex.htk.cpp: $(srcdir)/htk.l
	$(LEX) $(LFLAGS) $(AM_LFLAGS) --outfile=lex.htk.cpp $(srcdir)/htk.l 

