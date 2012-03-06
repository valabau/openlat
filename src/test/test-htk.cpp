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
#include <openlat/compat.h>
#include <openlat/rmarc.h>
#include <openlat/utils.h>
#include <openlat/htk-compiler.h>
#include <openlat/approx-shortest-distance.h>




using namespace fst;
using namespace openlat;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestHtk
#include <boost/test/unit_test.hpp>


struct Lattice {
};

BOOST_FIXTURE_TEST_SUITE(Htk, Lattice)

BOOST_AUTO_TEST_CASE(Htk)
{
  //LogLinearFst *fst = ReadHtk();
}

BOOST_AUTO_TEST_SUITE_END()
