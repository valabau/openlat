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
 * query.h
 *
 *  Created on: 20/02/2012
 *      Author: valabau
 */

#ifndef openlat_QUERY_H_
#define openlat_QUERY_H_

namespace openlat {

typedef fst::VectorFst<fst::LogArc> LogVectorFst;
typedef LogVectorFst::Arc::Label VLabel;

typedef struct _sSampleLabel{
  size_t sample;
  VLabel label;
  _sSampleLabel(size_t _sample = 0, VLabel _label = 0): sample(_sample), label(_label) {}
  bool operator==(const _sSampleLabel &other) const { return sample == other.sample and label == other.label; }
  bool operator<(const _sSampleLabel &other) const { return (sample != other.sample)?(sample < other.sample):(label < other.label); }
} SampleLabel;



typedef struct _sQuery {
  SampleLabel label;
  VLabel hyp;
  _sQuery(size_t _sample = 0, VLabel _label = 0, VLabel _hyp = fst::SymbolTable::kNoSymbol): label(_sample, _label), hyp(_hyp) {}
} Query;

// True if specified labels match (don't match) when M is
// true (false) following the T match_t type criteria
template <class A>
class QueryFilter {
public:
  typedef typename A::Label Label;

  QueryFilter(Label label, Label member): label_(label), member_(member) {};

  bool operator()(const A &arc) const {
    if (label_ != arc.ilabel) return true;
    else return member_ == arc.olabel;
  }

private:
  Label label_, member_;
};

typedef QueryFilter<fst::LogArc> LogQueryFilter;


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
    bool is_final = best->Final(state) != Arc::Weight::Zero();
    assert_bt((not is_final and best->NumArcs(state) == 1)
           or (is_final and best->NumArcs(state) == 0), "Unexpected number of arcs in 1-best fst\n");
    typename Arc::StateId nextstate = fst::kNoStateId;
    for (fst::ArcIterator<Fst> aiter(*best, state); !aiter.Done(); aiter.Next()) {
      typename Fst::Arc arc = aiter.Value();
      hyp.push_back(arc.olabel);
      scores.push_back(arc.weight.Value());
      state = arc.nextstate;
    }
    state = nextstate;
  }

  std::reverse(hyp.begin(), hyp.end());
  std::reverse(scores.begin(), scores.end());

  delete best;

  return score;
}

template<>
float ShortestPath<fst::LogArc>(const fst::Fst<fst::LogArc> &fst, vector<VLabel> &hyp, vector<float> &scores) {
  fst::MapFst<fst::LogArc, fst::StdArc, fst::LogToStdMapper> _fst(fst, fst::LogToStdMapper());
  return ShortestPath<fst::StdArc>(_fst, hyp, scores);
}



// WeightConvert class specialization that converts the weights.
// Mapper that changes nextstate to a dead state based on an ArcFilter
template <class A, class ArcFilter>
class FilterMapper {
 public:
  typedef A FromArc;
  typedef A ToArc;
  typedef ArcFilter Filter;
  typedef typename FromArc::Weight FromWeight;
  typedef typename ToArc::Weight ToWeight;

  FilterMapper(const Filter &filter, typename A::StateId dead_state): filter_(filter), dead_state_(dead_state) {};

  ToArc operator()(const FromArc &arc) const {
    return ToArc(arc.ilabel, arc.olabel,
                 arc.weight, (filter_(arc))?arc.nextstate:dead_state_);
  }

  fst::MapFinalAction FinalAction() const { return fst::MAP_NO_SUPERFINAL; }

  fst::MapSymbolsAction InputSymbolsAction() const { return fst::MAP_COPY_SYMBOLS; }

  fst::MapSymbolsAction OutputSymbolsAction() const { return fst::MAP_COPY_SYMBOLS;}

  uint64 Properties(uint64 props) const { return props; }

  void SetDeadState(typename A::StateId dead_state) { dead_state_ = dead_state; }

 private:
  const ArcFilter &filter_;
  typename A::StateId dead_state_;
};

typedef FilterMapper<fst::LogArc, LogQueryFilter> QueryMapper;

}

namespace std {
  namespace tr1 {
    template<> struct hash<openlat::SampleLabel> {
      size_t operator()(const openlat::SampleLabel& __x) const {
        return hash<size_t>()(__x.sample << 16) | hash<openlat::VLabel>()(__x.label);
      }
    };
  }
}


#endif /* openlat_QUERY_H_ */
