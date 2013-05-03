AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
AM_LDFLAGS = -L/usr/lib/fst -ldl -lfst -lfstscript #-lngram-fst 

bin_PROGRAMS = lat-compile lat-normalize lat-rmprefix lat-sumsuffixes lat-split-symbols lat-print-htk interactive-sequence-labeling #ngram-compile 
LDADD = ../lib/libopenlat.la -lfst -ldl $(BOOST_IOSTREAMS_LIB)

interactive_sequence_labeling_SOURCES = interactive-sequence-labeling.cpp
lat_compile_SOURCES = lat-compile.cpp
#ngram_compile_SOURCES = ngram-compile.cpp
lat_normalize_SOURCES = lat-normalize.cpp
lat_rmprefix_SOURCES = lat-rmprefix.cpp
lat_sumsuffixes_SOURCES = lat-sumsuffixes.cpp
lat_split_symbols_SOURCES = lat-split-symbols.cpp
lat_print_htk_SOURCES = lat-print-htk.cpp
