/*
 *   Copyright 2012, valabau
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 *
 * 
 * test-Htk.cpp
 *
 *  Created on: 17/02/2012
 *      Author: valabau
 */


#include <cmath>
#include <memory>

#include <fst/fstlib.h>
#include <openlat/config.h>
#include <openlat/compat.h>
#include <openlat/iofilter.h>




using namespace std;
using namespace fst;
using namespace openlat;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestHtk
#include <boost/test/unit_test.hpp>

string rw_word(string fn, string word) {
  ofilter ofl(fn.c_str());
  ofl << word;
  ofl.close();

  string res;
  ifilter ifl(fn.c_str());
  ifl >> res;

  return res;
}

struct FileFixture {
};

BOOST_FIXTURE_TEST_SUITE(File, FileFixture)

BOOST_AUTO_TEST_CASE(File)
{
  string res = rw_word("test-file.txt", "plain");
  cerr << "'" << res << "' == 'plain'\n";
  BOOST_CHECK(res == "plain");

#ifdef HAVE_BOOST_IOSTREAMS
  res = rw_word("test-file-gzipped.txt.gz", "gzipped");
  cerr << "'" << res << "' == 'gzipped'\n";
  BOOST_CHECK(res == "gzipped");

  res = rw_word("test-file-bzipped.txt.bz2", "bzipped");
  cerr << "'" << res << "' == 'bzipped'\n";
  BOOST_CHECK(res == "bzipped");
#endif
}

BOOST_AUTO_TEST_SUITE_END()
