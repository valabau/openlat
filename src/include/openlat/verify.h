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
 * verify.h
 *
 *  Created on: 06/03/2012
 *      Author: valabau
 */

#ifndef openlat_VERIFY_H_
#define openlat_VERIFY_H_

#include <fst/fst.h>

namespace openlat {

/** Verify that the FST is a Sequence Labeling FST and return the number of labels
 *
 *  A FST is a sequence labeling fst when all paths reaching a state have the same length.
 *  The number of labels, then, equals to the length of the paths to the final state.
 *
 *  @param fst a fst
 *  @return the length of the paths to the final state (number of labels) or (size_t)-1 if not a sequence labeling FST
 */
template<class Arc>
size_t VerifySequenceLabeling(const fst::Fst<Arc> &fst) {
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId StateId;

  if (not fst.Properties(fst::kTopSorted, true)) {
    LOG(ERROR) << "VerifySequenceLabeling: the fst must be topologically sorted to verify sequence labeling\n";
    return static_cast<size_t>(-1);
  }

  if (not fst.Properties(fst::kAcyclic, true)) {
    LOG(ERROR) << "VerifySequenceLabeling: the fst must be acyclic to verify sequence labeling\n";
    return static_cast<size_t>(-1);
  }

  // XXX: epsilons are treated as regular symbols
//  if (not fst.Properties(fst::kNoEpsilons, true)) {
//    LOG(ERROR) << "the fst must have no epsilons to verify sequence labeling\n";
//    return static_cast<size_t>(-1);
//  }

//  if (not fst.Properties(fst::kCoAccessible, true)) {
//    LOG(ERROR) << "the fst must be coaccessible to verify sequence labeling\n";
//    return static_cast<size_t>(-1);
//  }

  if (not Verify(fst)) return static_cast<size_t>(-1);

  // Count states
  StateId ns = 0;
  for (fst::StateIterator< fst::Fst<Arc> > siter(fst);
       !siter.Done();
       siter.Next())
    ++ns;

  vector<size_t> state_length(ns, static_cast<size_t>(-1));
  state_length[fst.Start()] = 0;

  size_t final_length = static_cast<size_t>(-1);

  for (fst::StateIterator< fst::Fst<Arc> > siter(fst);
       !siter.Done();
       siter.Next()) {
    StateId s = siter.Value();
    size_t this_length = state_length[s];
    if (this_length != static_cast<size_t>(-1)) {
      size_t na = 0;
      for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s);
           !aiter.Done();
           aiter.Next())
      {
        const Arc &arc = aiter.Value();
        if (arc.ilabel < 0) {
          LOG(ERROR) << "VerifySequenceLabeling: Fst input label ID of arc at position "
                     << na << " of state " << s << " is negative";
          return static_cast<size_t>(-1);
        } else if (arc.olabel < 0) {
          LOG(ERROR) << "VerifySequenceLabeling: Fst output label ID of arc at position "
                     << na << " of state " << s << " is negative";
          return static_cast<size_t>(-1);
        }

        if (state_length[arc.nextstate] == static_cast<size_t>(-1)) state_length[arc.nextstate] = this_length + 1;
        else if (state_length[arc.nextstate] != this_length + 1) {
          LOG(ERROR) << "VerifySequenceLabeling: Fst state " << state_length[arc.nextstate]
                     << " length from " << s << " is different from a previous state";
          return static_cast<size_t>(-1);
        }

        ++na;
      }
      if (!fst.Final(s).Member()) {
        LOG(ERROR) << "VerifySequenceLabeling: Fst final weight of state " << s << " is invalid";
        return static_cast<size_t>(-1);
      }
      if (fst.Final(s).Value() < Weight::Zero().Value()) {
        if (final_length == static_cast<size_t>(-1)) final_length = this_length;
        else if (this_length != final_length) {
          LOG(ERROR) << "VerifySequenceLabeling: Fst state " << s << " length is different from a previous final state length";
          return static_cast<size_t>(-1);
        }
      }
    }
  }

  return final_length;
}


}

#endif /* openlat_VERIFY_H_ */
