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
 * test-query.cpp
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
#include <openlat/query.h>
#include <openlat/approx-shortest-distance.h>




using namespace fst;
using namespace openlat;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestQuery
#include <boost/test/unit_test.hpp>


struct Lattice {
  LogVectorFst fst;


  Lattice() {
    BOOST_TEST_MESSAGE("initializing fst");

    fst.SetStart(fst.AddState());

    fst.AddArc(0, LogArc(0, 0, -log(0.4), 1));
    fst.AddArc(0, LogArc(0, 1, -log(0.6), 1));

    // Adds state 1 and its arc.
    fst.AddState();
    fst.AddArc(1, LogArc(1, 2, -log(0.5), 2));
    fst.AddArc(1, LogArc(1, 0, -log(0.5), 2));

    // Adds state 2 and set its final weight.
    fst.AddState();
    fst.SetFinal(2, -log(1.0));  // 1st arg is state ID, 2nd arg weigh

    Connect(&fst);
  }

  ~Lattice() {
      BOOST_TEST_MESSAGE("finalizing fst");
  }
};

typedef MapFst<LogArc, LogArc, QueryMapper> QueryFst;

BOOST_FIXTURE_TEST_SUITE(Query, Lattice)

BOOST_AUTO_TEST_CASE(removeArc)
{
  for (VLabel label = 0; label < 2; label++) {
    for (VLabel member = 0; member < 3; member++) {
      float distance, viterbi_distance;

      { // test arc remover
        LogVectorFst *fst1 = fst.Copy(true);
        BOOST_CHECK(Verify(*fst1));

        LogQueryFilter filter(label, member);
        RmArc<LogArc, LogQueryFilter>(fst1, RmArcOptions<LogArc, LogQueryFilter>(filter));

        for (StateIterator<LogVectorFst> siter(*fst1); !siter.Done(); siter.Next()) {
          for (ArcIterator<LogVectorFst> aiter(*fst1, siter.Value()); !aiter.Done(); aiter.Next()) {
            LogArc arc = aiter.Value();

            if (arc.ilabel == label) {
              BOOST_CHECK(arc.olabel == member);
            }
          }
        }

        distance = ShortestDistance(*fst1).Value();
        viterbi_distance = ApproxShortestDistance<LogArc, StdArc>(*fst1).Value();


        delete fst1;
      }

      {
        cerr << "distance:" << exp(-ShortestDistance(fst).Value()) << ":" << exp(-ApproxShortestDistance<LogArc, StdArc>(fst).Value()) << "\n";

        LogVectorFst::StateId dead_state = fst.AddState();

        LogQueryFilter filter(label, member);
        QueryMapper mapper(filter, dead_state);

        QueryFst fst1(fst, mapper);

        BOOST_CHECK(Verify(fst1));


        float dist = ShortestDistance(fst1).Value();
        float viterbi_dist = ApproxShortestDistance<LogArc, StdArc>(fst1).Value();
        BOOST_CHECK(distance == dist);
        BOOST_CHECK(viterbi_distance == viterbi_dist);
        cerr << "label:" << label << ", member:" << member << ", distance:" << exp(-distance) << ":" << exp(-dist)
             << ", vit:" << exp(-viterbi_distance) << ":" << exp(-viterbi_dist) << "\n";

        for (StateIterator<QueryFst> siter(fst1); !siter.Done(); siter.Next()) {
          for (ArcIterator<QueryFst> aiter(fst1, siter.Value()); !aiter.Done(); aiter.Next()) {
            LogArc arc = aiter.Value();

            if (arc.ilabel == label) {
              BOOST_CHECK(arc.olabel == member or arc.nextstate == dead_state);
            }
          }
        }

        fst.DeleteStates(vector<LogVectorFst::StateId>(1, dead_state));
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
