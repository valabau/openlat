AM_CPPFLAGS = $(WARNING_CXXFLAGS) $(COMPILANCE_CXXFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/src/include
AM_LDFLAGS =  $(BOOST_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

bin_PROGRAMS = test-query test-htk test-file test-sequence-labeling test-normalize
#LDADD = ../lib/libopenlat.la  
LDADD = ../lib/libopenlat.la -lfst -ldl -lz -lbz2  

test_query_SOURCES = test-query.cpp
test_htk_SOURCES = test-htk.cpp
test_file_SOURCES = test-file.cpp
test_sequence_labeling_SOURCES = test-sequence-labeling.cpp
test_normalize_SOURCES = test-normalize.cpp
