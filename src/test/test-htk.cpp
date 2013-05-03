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
#include <sstream>

#include <fst/fstlib.h>
#include <openlat/compat.h>
#include <openlat/rmarc.h>
#include <openlat/utils.h>
#include <openlat/htk-compiler.h>
#include <openlat/iofilter.h>
#include <openlat/approx-shortest-distance.h>




using namespace std;
using namespace fst;
using namespace openlat;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestHtk
#include <boost/test/unit_test.hpp>

extern int htk_debug;

struct LatticeTest {
  string file_content;

  LatticeTest() {
    file_content = "# Header (generated by SRILM)\n"
                    "VERSION=1.1\n"
                    "UTTERANCE=25730992H\n"
                    "NODES=173 LINKS=1513\n"
                    "# Nodes\n"
                    "I=0\n"
                    "I=1\n"
                    "# Links\n"
                    "J=0 S=0 E=1 W=Hello a=-10\n"
                    "J=1 S=0 E=1 W=Crocodile a=-12\n";
  }
};

BOOST_FIXTURE_TEST_SUITE(Htk, LatticeTest)

BOOST_AUTO_TEST_CASE(Htk)
{
  //LogLinearFst *fst = ReadHtk<LogLinearFst>();
  {
    istringstream is(file_content);
    MutableFst<LogArc> *fst = ReadHtkLogArc(is, "<string>");
    Verify(*fst);

    fst->Write("test-lat.fst");

    delete fst;
  }
  {
    ifilter is("src/test/dni.lat.bz2");
    MutableFst<LogArc> *fst = ReadHtkLogArc(is, "src/test/dni.lat.bz2");
    Verify(*fst);

    fst->Write("test-lat2.fst");

    delete fst;
  }
}

BOOST_AUTO_TEST_SUITE_END()
