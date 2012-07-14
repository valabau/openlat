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


namespace openlat {

template<class Arc>
float ShortestPath(const fst::Fst<Arc> &fst, vector<VLabel> &hyp, vector<float> &scores) {
  typedef fst::VectorFst<Arc> Fst;
  Fst *best = new Fst();

  ShortestPath(fst, best);

  float score = 0;
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
      nextstate = arc.nextstate;
    }
    state = nextstate;
  }

  delete best;

  return score;
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
    return __x->second.Value() < __y->second.Value();
  }
};


template <class Arc>
float ShortestHammingPath(const fst::Fst<Arc> &fst, vector<VLabel> &hyp, vector<float> &scores) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<Weight> forward;
  std::vector<Weight> backward;

  fst::ShortestDistance(fst, &forward);
  fst::ShortestDistance(fst, &backward, true);

  std::vector< std::map<Label, Weight> > probs;

  std::vector<size_t> position(NumStates(fst), static_cast<size_t>(-1));
  position[fst.Start()] = 0;

  // accumulate the posterior probabilities 
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();
    

    // ShortestDistance does not fill states that are not reachable at the end of the vector 
    if (s >= forward.size()) break; 
  

    if (position[s] != static_cast<size_t>(-1)) {
      const size_t pos = position[s];

      for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
        Arc arc = aiter.Value();

        if (position[arc.nextstate] == static_cast<size_t>(-1)) position[arc.nextstate] = pos + 1;
        else assert_bt(position[arc.nextstate] == pos + 1, "Unexpected state position. The fst seems not to be a sequence labeling problem");

        if (probs.size() < pos + 1) probs.resize(pos + 1);

        if (backward[arc.nextstate] != Weight::Zero() and arc.weight != Weight::Zero()) {
          typename std::map<Label, Weight>::iterator it = probs[pos].find(arc.olabel);
          if (it == probs.end()) {
            it = probs.insert(make_pair(arc.olabel, Weight::Zero())).first;
          }
          it->second = fst::Plus(it->second, fst::Divide(fst::Times(fst::Times(forward[s], arc.weight), backward[arc.nextstate]), backward[fst.Start()]) );
        }
      }
    }
  }

  hyp.clear();
  hyp.resize(probs.size());
  scores.clear();
  scores.resize(probs.size());

  Weight expected_num_errors = .0f;
  for (size_t i = 0; i < probs.size(); i++) {
    typename std::map<Label, Weight>::iterator it = std::max(probs[i].begin(), probs[i].end(), sort_by_score<Arc>());
    hyp[i] = it->first;
    scores[i] = - it->second.Value();
    expected_num_errors = fst::Plus( expected_num_errors, fst::Minus(Weight::One, it->second) );
  }

  return - expected_num_errors.Value();

}


}  // namespace openlat


#endif /* openlat_SEARCH_PATH_H_ */
