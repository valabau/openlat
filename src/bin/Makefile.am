AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
AM_LDFLAGS = -L/usr/lib/fst  

bin_PROGRAMS = lat-compile lat-normalize lat-rmprefix lat-sumsuffixes lat-split-symbols lat-print-htk lat-path-count lat-prune-arcs interactive-sequence-labeling #ngram-compile 
LDADD = ../lib/libopenlat.la -lfst -lfstscript $(BOOST_IOSTREAMS_LIB)  -ldl

interactive_sequence_labeling_SOURCES = interactive-sequence-labeling.cpp
lat_compile_SOURCES = lat-compile.cpp
#ngram_compile_SOURCES = ngram-compile.cpp
lat_normalize_SOURCES = lat-normalize.cpp
lat_rmprefix_SOURCES = lat-rmprefix.cpp
lat_sumsuffixes_SOURCES = lat-sumsuffixes.cpp
lat_split_symbols_SOURCES = lat-split-symbols.cpp
lat_print_htk_SOURCES = lat-print-htk.cpp
lat_path_count_SOURCES = lat-path-count.cpp
lat_prune_arcs_SOURCES = lat-prune-arcs.cpp
