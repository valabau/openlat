AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
AM_LDFLAGS = -lfst

bin_PROGRAMS = interactive-sequence-labeling lat-compile lat-normalize lat-rmprefix lat-sumsuffixes
LDADD = ../lib/libopenlat.la -lfst -ldl $(BOOST_IOSTREAMS_LIB)

interactive_sequence_labeling_SOURCES = interactive-sequence-labeling.cpp
lat_compile_SOURCES = lat-compile.cpp
lat_normalize_SOURCES = lat-normalize.cpp
lat_rmprefix_SOURCES = lat-rmprefix.cpp
lat_sumsuffixes_SOURCES = lat-sumsuffixes.cpp
