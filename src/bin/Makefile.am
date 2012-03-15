AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
AM_LDFLAGS = -lfst

bin_PROGRAMS = interactive-sequence-labeling lat-compile lat-normalize
LDADD = ../lib/libopenlat.la -lfst -ldl -lz -lbz2  

interactive_sequence_labeling_SOURCES = interactive-sequence-labeling.cpp
lat_compile_SOURCES = lat-compile.cpp
lat_normalize_SOURCES = lat-normalize.cpp
