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
 * approx-shortest-distance.h
 *
 *  Created on: 20/02/2012
 *      Author: valabau
 */

#ifndef openlat_APPROX_SHORTEST_DISTANCE_H_
#define openlat_APPROX_SHORTEST_DISTANCE_H_

#include <fst/shortest-distance.h>

namespace openlat {

// Approximate versions of the shortest-distance algorithms.
// It works by converts LogWeights to TropicalWeights.
// For information on each function, see the original shortest-distance documentation.

template<class FromArc, class ToArc, class Queue, class ArcFilter>
void ApproxShortestDistance(
    const Fst<FromArc> &fst,
    vector<typename ToArc::Weight> *distance,
    const ShortestDistanceOptions<FromArc, Queue, ArcFilter> &opts) {
  typedef WeightConvertMapper<FromArc, ToArc> WeightMapper;
  typedef MapFst<FromArc, ToArc, WeightMapper> ApproxFst;
  ShortestDistance<ToArc, Queue, ArcFilter>(ApproxFst(fst, WeightMapper()), distance, opts);
}

template <class FromArc, class ToArc>
void ApproxShortestDistance(const Fst<FromArc> &fst,
                      vector<typename ToArc::Weight> *distance,
                      bool reverse = false,
                      float delta = kDelta) {
  typedef WeightConvertMapper<FromArc, ToArc> WeightMapper;
  typedef MapFst<FromArc, ToArc, WeightMapper> ApproxFst;
  ShortestDistance<ToArc>(ApproxFst(fst, WeightMapper()), distance, reverse, delta);
}


// Return the sum of the weight of all successful paths in an FST, i.e.,
// the shortest-distance from the initial state to the final states.
template <class FromArc, class ToArc>
typename ToArc::Weight ApproxShortestDistance(const Fst<FromArc> &fst, float delta = kDelta) {
  typedef WeightConvertMapper<FromArc, ToArc> WeightMapper;
  typedef MapFst<FromArc, ToArc, WeightMapper> ApproxFst;
  return ShortestDistance<ToArc>(ApproxFst(fst, WeightMapper()), delta);
}


}  // namespace openfst



#endif /* openlat_APPROX_SHORTEST_DISTANCE_H_ */
