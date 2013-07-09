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

#ifndef openlat_PATH_COUNT_H_
#define openlat_PATH_COUNT_H_

#include <openlat/normalize.h>

namespace openlat {

template<class Arc>
void StatePathCount(const fst::Fst<Arc> &fst, std::vector<fst::LogArc::Weight> *distance) {
  typedef fst::VectorFst<fst::LogArc> Fst;
  typedef typename Arc::Weight Weight;
  Fst *pfst = new Fst();

  fst::ArcMap(fst, pfst, fst::WeightConvertMapper<Arc, Fst::Arc>());
  fst::ArcMap(pfst, fst::RmWeightMapper<Fst::Arc>());

  
  ShortestDistance(*pfst, distance, true);

  delete pfst;
}


template<class Arc>
float PathCount(const fst::Fst<Arc> &fst) {
  vector<fst::LogArc::Weight> distance;
  StatePathCount(fst, &distance);
  return exp(-distance[fst.Start()].Value());
}


}  // namespace openlat


#endif /* openlat_PATH_COUNT_H_ */
