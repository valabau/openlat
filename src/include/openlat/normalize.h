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
 * normalize.h
 *
 *  Created on: 15/03/2012
 *      Author: valabau
 */

#ifndef openlat_NORMALIZE_H_
#define openlat_NORMALIZE_H_

#include <fst/shortest-distance.h>

namespace openlat {

// Mapper to (right) multiply a constant to all weights.
template <class A>
struct PowerMapper {
  typedef A FromArc;
  typedef A ToArc;
  typedef typename A::Weight Weight;

  explicit PowerMapper(float w) : weight_(w) {}

  A operator()(const A &arc) const {
    if (arc.weight == Weight::Zero())
      return arc;
    Weight w = Weight(arc.weight.Value() * weight_);
    return A(arc.ilabel, arc.olabel, w, arc.nextstate);
  }

  fst::MapFinalAction FinalAction() const { return fst::MAP_NO_SUPERFINAL; }

  fst::MapSymbolsAction InputSymbolsAction() const { return fst::MAP_COPY_SYMBOLS; }

  fst::MapSymbolsAction OutputSymbolsAction() const { return fst::MAP_COPY_SYMBOLS;}

  uint64 Properties(uint64 props) const {
    return props & fst::kWeightInvariantProperties;
  }

 private:
  float weight_;
};


template<class Arc>
bool VerifyProbabilistic(const fst::Fst<Arc> &fst, float delta = fst::kDelta) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();
    Weight sum = Weight::Zero();

    size_t n = 0;
    for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      sum = fst::Plus(sum, aiter.Value().weight);
      n++;
    }
    sum = fst::Plus(sum, fst.Final(s));

    if (delta != 0) {
      if (not fst::ApproxEqual(sum, Weight::One(), delta)) return false;
    }
    else {
      if (sum != Weight::One()) return false;
    }
  }

  return true;
}

template<class Arc>
void Normalize(fst::MutableFst<Arc> *fst) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<Weight> backward;
  fst::ShortestDistance(*fst, &backward, true);

  for (fst::StateIterator<fst::MutableFst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    // normalize arcs
    for (fst::MutableArcIterator < fst::MutableFst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      arc.weight = fst::Divide(fst::Times(arc.weight, backward[arc.nextstate]), backward[s]);
      aiter.SetValue(arc);
    }

    // normalize final probability
    if (fst->Final(s) != Weight::Zero()) {
      fst->SetFinal(s, fst::Divide(fst->Final(s), backward[s]));
    }
  }

}

template<class Arc>
void DeterminizeAndNormalize(const fst::Fst<Arc> &ifst,
                             fst::MutableFst<Arc> *ofst,
                             const fst::DeterminizeOptions<Arc> &opts
                             = fst::DeterminizeOptions<Arc>())
{
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  fst::VectorFst<Arc> fst(ifst);

  if (fst.Properties(fst::kEpsilons, true) & fst::kEpsilons) {
    fst::RmEpsilon(&fst);
  }

  fst::EncodeMapper<Arc> encode_mapper(fst::kEncodeLabels, fst::ENCODE);
  fst::Encode(&fst, &encode_mapper);

  fst::Determinize(fst, ofst, opts);
  fst::Minimize(ofst);

  fst::Decode(ofst, encode_mapper);
  Normalize(ofst);
}


template <class W>
inline float Entropy(const W &w, const W &fwd) {
  if (not w.Member() or w == W::Zero()) return w.Value();
  return -w.Value() * exp(-w.Value() -fwd.Value());
}


template<class Arc>
float Entropy(const fst::Fst<Arc> &fst) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<Weight> forward;
  std::vector<Weight> backward;

  fst::ShortestDistance(fst, &forward);
  fst::ShortestDistance(fst, &backward, true);

  float perplexity = .0;
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      perplexity += Entropy(arc.weight, forward[s]);
    }

    // add entropy for final state
    if (fst.Final(s) != Weight::Zero()) {
      perplexity += Entropy(fst.Final(s), forward[s]);
    }

  }

  return perplexity;
}


}

#endif /* openlat_NORMALIZE_H_ */
