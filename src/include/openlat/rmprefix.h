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
 * rmprefix.h
 *
 *  Created on: 12/04/2012
 *      Author: valabau
 */

#ifndef openlat_RMPREFIX_H_
#define openlat_RMPREFIX_H_

#include <fst/mutable-fst.h>
#include <tr1/unordered_map>

namespace openlat {

template<class Arc>
bool RmPrefix(fst::MutableFst<Arc> *fst, size_t prefix_length) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  if (prefix_length == 0) return true;

  Weight push_weight = Weight::One();
  std::vector<StateId> states;
  states.reserve(prefix_length);
  states.push_back(fst->Start());

  for (size_t p = 0; p < prefix_length; p++) {
    size_t narcs = 0;
    for (fst::ArcIterator < fst::Fst<Arc> > aiter(*fst, states.back()); !aiter.Done(); aiter.Next()) {
      const Arc &arc = aiter.Value();
      push_weight = Times(push_weight, arc.weight);
      states.push_back(arc.nextstate);
      narcs++;
    }

    if (narcs != 1) {
      return false;
    }
  }

  fst->SetStart(states.back());

  for (fst::MutableArcIterator < fst::MutableFst<Arc> > aiter(fst, states.back()); !aiter.Done(); aiter.Next()) {
    Arc arc = aiter.Value();
    arc.weight = Times(push_weight, arc.weight);
    aiter.SetValue(arc);
  }
  if (fst->Final(states.back()) != Weight::Zero()) {
    fst->SetFinal(states.back(), Times(push_weight, fst->Final(states.back())));
  }

  states.pop_back();
  fst->DeleteStates(states);

  return true;
}


template<class Arc>
void SumSuffixes(const fst::Fst<Arc> &fst, fst::MutableFst<Arc> *wordlist) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  typedef std::tr1::unordered_map<typename Arc::Label, typename Arc::Weight> LabelWeightMap;

  std::vector<Weight> backward;
  fst::ShortestDistance(fst, &backward, true);

  Weight norm = backward[fst.Start()];

  LabelWeightMap wordweights;

  for (fst::ArcIterator < fst::Fst<Arc> > aiter(fst, fst.Start()); !aiter.Done(); aiter.Next()) {
    const Arc &arc = aiter.Value();

    typename LabelWeightMap::iterator weight_it = wordweights.find(arc.ilabel);
    if (weight_it != wordweights.end()) {
      wordweights[arc.ilabel] = Plus(wordweights[arc.ilabel], Divide(arc.weight, norm));
    }
    else {
      wordweights.insert(std::make_pair(arc.ilabel, Divide(arc.weight, norm)));
    }
  }
  if (fst.Final(fst.Start()) != Weight::Zero()) {
    wordweights.insert(std::make_pair(0, Divide(fst.Final(fst.Start()), norm)));
  }

  typename LabelWeightMap::iterator w_it;
  StateId start = wordlist->AddState();
  wordlist->SetStart(start);

  StateId end = wordlist->AddState();
  wordlist->SetFinal(end, Weight::One());
  for (w_it = wordweights.begin(); w_it != wordweights.end(); ++w_it) {
    wordlist->AddArc(start, Arc(w_it->first, w_it->first, w_it->second, end));
  }

  wordlist->SetInputSymbols(fst.InputSymbols());
  wordlist->SetOutputSymbols(fst.OutputSymbols());
}

}

#endif /* openlat_RMPREFIX_H_ */
