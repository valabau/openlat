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
 * search-path.h
 *
 *  Created on: 20/02/2012
 *      Author: valabau
 */

#ifndef openlat_SEARCH_PATH_H_
#define openlat_SEARCH_PATH_H_

#include <openlat/normalize.h>

namespace openlat {

template<class Arc>
float ShortestPath(const fst::Fst<Arc> &fst, vector<VLabel> &hyp, vector<float> &scores) {
  typedef fst::VectorFst<Arc> Fst;
  typedef typename Arc::Weight Weight;
  Fst *best = new Fst();

  ShortestPath(fst, best);

  Weight score = Weight::One();
  hyp.clear();
  scores.clear();

  typename Arc::StateId state = best->Start();
  while (state != fst::kNoStateId) {
    // cerr << "state: " << state << "\n";
    bool is_final = best->Final(state) != Arc::Weight::Zero();
    assert_bt((not is_final and best->NumArcs(state) == 1)
           or (is_final and best->NumArcs(state) == 0), "Unexpected number of arcs in 1-best fst\n");
    typename Arc::StateId nextstate = fst::kNoStateId;
    for (fst::ArcIterator<Fst> aiter(*best, state); !aiter.Done(); aiter.Next()) {
      typename Fst::Arc arc = aiter.Value();
      // cerr << "arc: " << arc.ilabel << " " << arc.olabel << " " << arc.weight.Value()  << " " << arc.nextstate << "\n";
      hyp.push_back(arc.olabel);
      scores.push_back(arc.weight.Value());
      score = fst::Times(score, arc.weight);
      nextstate = arc.nextstate;
    }
    state = nextstate;
  }

  delete best;


  return -score.Value();
}

template<>
float ShortestPath<fst::LogArc>(const fst::Fst<fst::LogArc> &fst, vector<VLabel> &hyp, vector<float> &scores) {
  fst::MapFst<fst::LogArc, fst::StdArc, fst::LogToStdMapper> _fst(fst, fst::LogToStdMapper());
  return ShortestPath<fst::StdArc>(_fst, hyp, scores);
}

/**
 * structure function to sort probs
 */
template <class Arc>
struct sort_by_score {
  size_t operator()(const std::pair<typename Arc::Label, typename Arc::Weight>& __x, const std::pair<typename Arc::Label, typename Arc::Weight>& __y) const {
    return __x.second.Value() > __y.second.Value();
  }
};


/**
 * Wrapper to use shortest path as template arguments
 */
template<class Arc>
struct ShortestPathFunc {
  float operator()(const fst::Fst<fst::LogArc> &fst, vector<VLabel> &hyp, vector<float> &scores) const {
    return ShortestPath<Arc>(fst, hyp, scores);
  }
};


template <class Arc>
void ExpectedHammingProbs(const fst::Fst<Arc> &fst, std::vector< std::map<typename Arc::Label, typename Arc::Weight> >& probs) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<Weight> forward;
  std::vector<Weight> backward;

  fst::ShortestDistance(fst, &forward);
  fst::ShortestDistance(fst, &backward, true);

  std::vector<size_t> position(NumStates(fst), static_cast<size_t>(-1));
  position[fst.Start()] = 0;

  // accumulate the posterior probabilities
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();


    // ShortestDistance does not fill states that are not reachable at the end of the vector
    if (static_cast<size_t>(s) >= forward.size()) break;


    if (position[s] != static_cast<size_t>(-1)) {
      const size_t pos = position[s];

      for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
        Arc arc = aiter.Value();

        if (position[arc.nextstate] == static_cast<size_t>(-1)) position[arc.nextstate] = pos + 1;
        else assert_bt(position[arc.nextstate] == pos + 1, "Unexpected state position. The fst seems not to be a sequence labeling problem");

        if (probs.size() < pos + 1) probs.resize(pos + 1);

        if (backward[arc.nextstate] != Weight::Zero() and arc.weight != Weight::Zero()) {
          typename std::map<Label, Weight>::iterator it = probs[pos].find(arc.olabel);
          if (it == probs[pos].end()) {
            it = probs[pos].insert(make_pair(arc.olabel, Weight::Zero())).first;
          }
          it->second = fst::Plus(it->second, fst::Divide(fst::Times(fst::Times(forward[s], arc.weight), backward[arc.nextstate]), backward[fst.Start()]) );
        }
      }
    }
  }

}



template <class Arc>
float ShortestHammingPath(const fst::Fst<Arc> &fst, vector<VLabel> &hyp, vector<float> &scores) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector< std::map<Label, Weight> > probs;
  ExpectedHammingProbs(fst, probs);

  hyp.clear();
  hyp.resize(probs.size());
  scores.clear();
  scores.resize(probs.size());

  Weight expected_num_ok = Weight::Zero();
  for (size_t i = 0; i < probs.size(); i++) {
    typename std::map<Label, Weight>::iterator it = std::max_element(probs[i].begin(), probs[i].end(), sort_by_score<Arc>());
    hyp[i] = it->first;
    scores[i] = - it->second.Value();
    expected_num_ok = fst::Plus(expected_num_ok, it->second);
  }

  //XXX: convert from expected accuracy to expected errors
  return to_float(expected_num_ok);

}

/**
 * Wrapper to use shortest path as template arguments
 */
template<class Arc>
struct ShortestHammingPathFunc {
  float operator()(const fst::Fst<fst::LogArc> &fst, vector<VLabel> &hyp, vector<float> &scores) const {
    return ShortestPath<Arc>(fst, hyp, scores);
  }
};



}  // namespace openlat


#endif /* openlat_SEARCH_PATH_H_ */
