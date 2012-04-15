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
 * test-sequence-labeling.cpp
 *
 *  Created on: 06/03/2012
 *      Author: valabau
 */




#include <cmath>
#include <memory>

#include <fst/fstlib.h>
#include <openlat/compat.h>
#include <openlat/rmarc.h>
#include <openlat/utils.h>
#include <openlat/query.h>
#include <openlat/verify.h>
#include <openlat/interactive-sequence-labeling.h>




using namespace fst;
using namespace openlat;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestSequenceLabeling
#include <boost/test/unit_test.hpp>

struct QueryResult {
  VLabel label;
  VLabel member;
  vector<VLabel> result;
};


struct Lattice {
  LogVectorFst fst;
  vector<VLabel> ref;
  vector<QueryResult> queries;


  Lattice() {
    BOOST_TEST_MESSAGE("initializing fst");
    srand(time(NULL));

    fst.SetStart(fst.AddState());

    fst.AddArc(0, LogArc(0, 0, -log(0.4), 1));
    fst.AddArc(0, LogArc(0, 1, -log(0.6), 1));

    // Adds state 1 and its arc.
    fst.AddState();
    fst.AddArc(1, LogArc(1, 2, -log(0.6), 2));
    fst.AddArc(1, LogArc(1, 0, -log(0.4), 2));

    // Adds state 2 and set its final weight.
    fst.AddState();
    fst.SetFinal(2, -log(1.0));  // 1st arg is state ID, 2nd arg weigh

    ref.push_back(0);
    ref.push_back(2);

    Connect(&fst);

    QueryResult q;
    q.result.resize(2);

    q.label = 0; q.member = 0; q.result[0] = 0; q.result[1] = 2;
    queries.push_back(q);

    q.label = 0; q.member = 1; q.result[0] = 1; q.result[1] = 2;
    queries.push_back(q);

    q.label = 1; q.member = 2; q.result[0] = 1; q.result[1] = 2;
    queries.push_back(q);

    q.label = 1; q.member = 0; q.result[0] = 1; q.result[1] = 0;
    queries.push_back(q);

  }

  ~Lattice() {
      BOOST_TEST_MESSAGE("finalizing fst");
  }
};

BOOST_FIXTURE_TEST_SUITE(SequenceLabeling, Lattice)

BOOST_AUTO_TEST_CASE(verifyFsts)
{
  BOOST_CHECK(Verify(fst));

  size_t len = VerifySequenceLabeling(fst);
  BOOST_CHECK(len == 2);

  { // state with different path lengths
    LogVectorFst *bad_fst = fst.Copy(true);
    bad_fst->AddArc(0, LogArc(0, 1, -log(0.6), 2));

    len = VerifySequenceLabeling(*bad_fst);
    BOOST_CHECK(len == static_cast<size_t>(-1));
    delete bad_fst;
  }

  { // fst with cycles
    LogVectorFst *bad_fst = fst.Copy(true);
    bad_fst->AddArc(1, LogArc(0, 1, -log(0.6), 1));

    len = VerifySequenceLabeling(*bad_fst);
    BOOST_CHECK(len == static_cast<size_t>(-1));
    delete bad_fst;
  }
}

BOOST_AUTO_TEST_CASE(mapper)
{
  BOOST_CHECK(Verify(fst));

  // obtain the shortest path with the restriction (label, member)
  vector<VLabel> lhyp;
  vector<float> lscores;
  ShortestPath(fst, lhyp, lscores);

  // cerr << syms_to_string(lhyp, 0) << "\n";
  BOOST_CHECK(lhyp.size() == 2);
  BOOST_CHECK(lhyp[0] == 1);
  BOOST_CHECK(lhyp[1] == 2);

  LogVectorFst::StateId dead_state = fst.AddState();
  for (vector<QueryResult>::const_iterator it = queries.begin(); it < queries.end(); ++it) {
    typedef LogArc Arc;
    typedef LogQueryFilter Filter;
    typedef FilterMapper<Arc, Filter> Mapper;

    const QueryResult &q = *it;
    Filter filter(q.label, q.member);
    Mapper qmapper(filter, dead_state);
    MapFst<Arc, Arc, Mapper> _fst(fst, qmapper);

    // obtain the shortest path with the restriction (label, member)
    lhyp.clear();
    lscores.clear();
    ShortestPath(_fst, lhyp, lscores);

    BOOST_CHECK(lhyp.size() == q.result.size());
    for (size_t i = 0; i < q.result.size(); i++) {
      BOOST_CHECK(lhyp[i] == q.result[i]);
    }
  }
  fst.DeleteStates(vector<LogVectorFst::StateId>(1, dead_state));

}

BOOST_AUTO_TEST_CASE(oracle)
{
  BOOST_CHECK(Verify(fst));

  vector<VectorFst<LogArc> *> fsts(1, &fst);
  vector<vector<VLabel> > refs(1, ref);

  LocalSystem<LogArc, LogQueryFilter, RecomputeExpectation<LogArc, LogQueryFilter> > system(fsts);
  Oracle oracle(refs);

  oracle.evaluate(&system);

}

BOOST_AUTO_TEST_SUITE_END()
