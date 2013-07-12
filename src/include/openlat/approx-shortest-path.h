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
 * approx-shortest-path.h
 *
 *  Created on: 20/02/2012
 *      Author: valabau
 */

#ifndef openlat_APPROX_SHORTEST_PATH_H_
#define openlat_APPROX_SHORTEST_PATH_H_

#include <fst/fst.h>

namespace openlat {

// Approximate versions of the shortest-path algorithms.
// It works by converts LogWeights to TropicalWeights.
// For information on each function, see the original shortest-path documentation.


template<class FromArc, class ToArc, class Queue, class ArcFilter>
void ApproxShortestPath(const fst::Fst<FromArc> &ifst, fst::MutableFst<ToArc> *ofst,
                  vector<typename ToArc::Weight> *distance,
                  fst::ShortestPathOptions<ToArc, Queue, ArcFilter> &opts) {
  typedef fst::WeightConvertMapper<FromArc, ToArc> WeightMapper;
  typedef fst::MapFst<FromArc, ToArc, WeightMapper> ApproxFst;
  fst::ShortestPath(ApproxFst(fst, WeightMapper), ofst, distance, opts);
}


template<class FromArc, ToArc>
void ApproxShortestPath(const fst::Fst<FromArc> &ifst, fst::MutableFst<ToArc> *ofst,
                  size_t n = 1, bool unique = false,
                  bool first_path = false,
                  typename ToArc::Weight weight_threshold = Arc::Weight::Zero(),
                  typename ToArc::StateId state_threshold = fst::kNoStateId) {
  typedef fst::WeightConvertMapper<FromArc, ToArc> WeightMapper;
  typedef fst::MapFst<FromArc, ToArc, WeightMapper> ApproxFst;
  fst::ShortestPath(ApproxFst(fst, WeightMapper), ofst, n, unique, first_path, weight_threshold, state_threshold);
}

}  // namespace openlat


#endif /* openlat_APPROX_SHORTEST_PATH_H_ */
