AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) -I$(top_srcdir)/src/include
AM_LDFLAGS = -L/usr/lib/fst -Wl,--no-as-needed 

bin_PROGRAMS = lat-compile lat-normalize lat-posteriors lat-rmprefix lat-sumsuffixes lat-split-symbols lat-print-htk lat-print-thot lat-path-count lat-prune-arcs lat-set-length interactive-sequence-labeling lat-confusion-network # ngram-compile 
LDADD =  -ldl ../lib/libopenlat.la -lfst -lfstscript $(BOOST_IOSTREAMS_LIB)

interactive_sequence_labeling_SOURCES = interactive-sequence-labeling.cpp
lat_compile_SOURCES = lat-compile.cpp
#ngram_compile_SOURCES = ngram-compile.cpp
lat_normalize_SOURCES = lat-normalize.cpp
lat_posteriors_SOURCES = lat-posteriors.cpp
lat_rmprefix_SOURCES = lat-rmprefix.cpp
lat_sumsuffixes_SOURCES = lat-sumsuffixes.cpp
lat_split_symbols_SOURCES = lat-split-symbols.cpp
lat_print_htk_SOURCES = lat-print-htk.cpp
lat_print_thot_SOURCES = lat-print-thot.cpp
lat_path_count_SOURCES = lat-path-count.cpp
lat_prune_arcs_SOURCES = lat-prune-arcs.cpp
lat_set_length_SOURCES = lat-set-length.cpp
lat_confusion_network_SOURCES = lat-confusion-network.cpp
