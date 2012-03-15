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

template<class Arc>
void Normalize(fst::MutableFst<Arc> *fst) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::vector<Weight> backward;

  fst::ShortestDistance(*fst, &backward, true);

  size_t n_arcs = 0;
  for (fst::StateIterator<fst::MutableFst<Arc> > siter(*fst); !siter.Done(); siter.Next()) n_arcs += fst->NumArcs(siter.Value());

  for (fst::StateIterator<fst::MutableFst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    for (fst::MutableArcIterator < fst::MutableFst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      arc.weight = fst::Divide(fst::Times(arc.weight, backward[arc.nextstate]), backward[s]);
      aiter.SetValue(arc);
    }
  }

}


}

#endif /* openlat_NORMALIZE_H_ */
