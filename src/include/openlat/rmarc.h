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

namespace openlat {

#include <fst/mutable-fst.h>
using namespace fst;


template <class Arc, class ArcFilter>
struct RmArcOptions {
  ArcFilter arc_filter;  // Arc filter (e.g., limit to only epsilon graph)
  RmArcOptions(ArcFilter filt): arc_filter(filt) {}
};


template<class Arc, class ArcFilter>
void RmArc(MutableFst<Arc> *fst, const RmArcOptions<Arc, ArcFilter> &opts) {
  typedef typename Arc::StateId StateId;

  StateId dead = fst->AddState();
  for (StateIterator<MutableFst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    for (MutableArcIterator < MutableFst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      if (not opts.arc_filter(arc)) {
        arc.nextstate = dead;
        aiter.SetValue(arc);
      }
    }
  }

  fst->DeleteStates(vector<StateId>(1, dead));
}


}

#endif /* openlat_RMARC_H_ */
