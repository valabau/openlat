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
 * rmarcs.h
 *
 *  Created on: 17/02/2012
 *      Author: valabau
 */

#ifndef openlat_RMARC_H_
#define openlat_RMARC_H_

#include <fst/mutable-fst.h>

namespace openlat {

template <class Arc, class ArcFilter>
struct RmArcOptions {
  ArcFilter arc_filter;  // Arc filter (e.g., limit to only epsilon graph)
  RmArcOptions(ArcFilter filt): arc_filter(filt) {}
};


template<class Arc, class ArcFilter>
void RmArc(fst::MutableFst<Arc> *fst, const RmArcOptions<Arc, ArcFilter> &opts) {
  typedef typename Arc::StateId StateId;

  StateId dead = fst->AddState();
  for (fst::StateIterator<fst::MutableFst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    for (fst::MutableArcIterator < fst::MutableFst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      if (not opts.arc_filter(arc)) {
        arc.nextstate = dead;
        aiter.SetValue(arc);
      }
    }
  }

  // fst->DeleteStates(vector<StateId>(1, dead));
  fst::Connect(fst);
}

void PruneLogArc(fst::MutableFst<fst::LogArc> *ifst, fst::StdArc::Weight threshold) {
  fst::VectorFst<fst::StdArc> stdfst;
  fst::ArcMap(*ifst, &stdfst, fst::WeightConvertMapper<fst::LogArc, fst::StdArc>());
  fst::Prune(&stdfst, threshold);
  fst::ArcMap(stdfst, ifst, fst::WeightConvertMapper<fst::StdArc, fst::LogArc>());
}

template<class Arc>
void ArcPosteriors(const fst::Fst<Arc> &fst, std::vector<typename Arc::Weight> *posteriors) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<Weight> forward;
  fst::ShortestDistance(fst, &forward);

  std::vector<Weight> backward;
  fst::ShortestDistance(fst, &backward, true);

  posteriors->clear();
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    // normalize arcs
    for (fst::ArcIterator < fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      posteriors->push_back(fst::Times(forward[s], fst::Times(arc.weight, backward[arc.nextstate])));
    }
  }
}

template<class Arc>
void PruneArcs(fst::MutableFst<Arc> *fst, typename Arc::Weight threshold, typename Arc::Weight zero = Arc::Weight::Zero()) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<typename Arc::Weight> posteriors;
  ArcPosteriors(*fst, &posteriors);
  

  Weight max_post = zero;
  for (size_t i = 0; i < posteriors.size(); ++i) {
    if (posteriors[i].Value() < max_post.Value()) max_post = posteriors[i];
  }

  threshold = fst::Times(threshold, max_post);
  //cerr << "max = " << max_post.Value() << "\n";
  //cerr << "thr = " << threshold.Value() << "\n";

  typename std::vector<Weight>::const_iterator post_it = posteriors.begin();
  for (fst::StateIterator<fst::MutableFst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    // set arc values to zero
    for (fst::MutableArcIterator < fst::MutableFst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next(), ++post_it) {
      Weight post = *post_it;
      //cerr << "post = " << post.Value() << "\n";
      if (post.Value() > threshold.Value()) {
        Arc arc = aiter.Value();
        arc.weight = zero; 
        aiter.SetValue(arc);
      }
    }
  }
  fst::Connect(fst);
}




}

#endif /* openlat_RMARC_H_ */
