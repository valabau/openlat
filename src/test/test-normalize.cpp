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
 * test-normalize.cpp
 *
 *  Created on: 15/03/2012
 *      Author: valabau
 */



#include <cmath>
#include <memory>

#include <fst/fstlib.h>
#include <openlat/compat.h>
#include <openlat/rmarc.h>
#include <openlat/utils.h>
#include <openlat/normalize.h>




using namespace fst;
using namespace openlat;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestNormalize
#include <boost/test/unit_test.hpp>

typedef fst::VectorFst<fst::LogArc> LogVectorFst;


struct Lattice {
  LogVectorFst fst;
  vector<LogWeight> weights;

  Lattice() {
    BOOST_TEST_MESSAGE("initializing fst");



    fst.SetStart(fst.AddState());

    fst.AddArc(0, LogArc(0, 0, -log(2*0.4), 1));
    fst.AddArc(0, LogArc(0, 1, -log(2*0.6), 1));

    // Adds state 1 and its arc.
    fst.AddState();
    fst.AddArc(1, LogArc(1, 2, -log(2*0.5), 2));
    fst.AddArc(1, LogArc(1, 0, -log(2*0.5), 2));

    // Adds state 2 and set its final weight.
    fst.AddState();
    fst.SetFinal(2, -log(1.0));  // 1st arg is state ID, 2nd arg weigh

    weights.push_back(-log(0.4));
    weights.push_back(-log(0.6));
    weights.push_back(-log(0.5));
    weights.push_back(-log(0.5));

    Connect(&fst);
  }

  ~Lattice() {
      BOOST_TEST_MESSAGE("finalizing fst");
  }
};

BOOST_FIXTURE_TEST_SUITE(NormalizeTest, Lattice)

BOOST_AUTO_TEST_CASE(normalize)
{
  Verify(fst);
  fst.Write("unormalized.fst");

  Normalize(&fst);
  float entropy = Entropy(fst);
  cerr << "entropy: " << entropy << "\n";
  cerr << "perplexity: " << exp(entropy) << "\n";

  vector<LogWeight>::const_iterator wit = weights.begin();

  for (StateIterator<MutableFst<LogArc> > siter(fst); !siter.Done(); siter.Next()) {
    LogArc::StateId s = siter.Value();

    for (MutableArcIterator<MutableFst<LogArc> > aiter(&fst, s); !aiter.Done(); aiter.Next()) {
      LogArc arc = aiter.Value();
      BOOST_CHECK(arc.weight == *wit);
      ++wit;
    }
  }
  fst.Write("normalized.fst");

  BOOST_CHECK(VerifyProbabilistic(fst, 0));
}

BOOST_AUTO_TEST_SUITE_END()
