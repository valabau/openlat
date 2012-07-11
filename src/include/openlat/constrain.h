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
 * constrain.h
 *
 *  Created on: 20/02/2012
 *      Author: valabau
 */

#ifndef openlat_CONSTRAIN_H_
#define openlat_CONSTRAIN_H_

#include <fst/fst.h>

namespace openlat {

template <typename Label>
class FstConstraint {
public:
  FstConstraint(): match_symbol_(fst::kNoLabel) {}
  ~FstConstraint() {}
  bool isAllowed(const Label& symbol) const {
    if (match_symbol_ != fst::kNoLabel) return symbol == match_symbol_;
    else return exclude_symbols_.find(symbol) == exclude_symbols_.end();
  }
  void clear() {
    match_symbol_ = fst::kNoLabel;
    exclude_symbols_.clear();
  }
  void setMatch(const Label& symbol) {
    match_symbol_ = symbol;
  }
  Label getMatch() const {
    return match_symbol_;
  }
  void addExclude(const Label& symbol) {
    exclude_symbols_.insert(symbol);
  }

private:
  Label match_symbol_;
  std::set<Label> exclude_symbols_;
};



// True if the constraints are not broken
// true (false) following the T match_t type criteria
template <class A>
class ConstraintFilter {
public:
  typedef typename A::Label Label;
  typedef typename std::map<Label, FstConstraint<Label> > ConstraintMap;

  ConstraintFilter() {};

  bool operator()(const A &arc) const {
    typename ConstraintMap::const_iterator it = constraints_.find(arc.ilabel);
    if (it == constraints_.end()) return true;
    else return it->second.isAllowed(arc.olabel);
  }
  void clear() {
    typename ConstraintMap::const_iterator it = constraints_.begin();
    while (it != constraints_.end()) {
      it->second.clear();
      ++it;
    }
  }
  void setMatch(const Label& label, const Label& symbol) {
    Label _label = label + 1; // XXX: sum 1 since label indexes start at 1
    typename ConstraintMap::iterator it = constraints_.find(_label);
    if (it == constraints_.end()) {
      it = constraints_.insert(std::make_pair(_label, FstConstraint<Label>())).first;
    }
    it->second.setMatch(symbol);
  }

  Label getMatch(const Label& label) const {
    Label _label = label + 1; // XXX: sum 1 since label indexes start at 1
    typename ConstraintMap::const_iterator it = constraints_.find(_label);
    if (it == constraints_.end()) return fst::kNoLabel;
    else return it->second.getMatch();
  }

  void addExclude(const Label& label, const Label& symbol) {
    Label _label = label + 1; // XXX: sum 1 since label indexes start at 1
    typename ConstraintMap::iterator it = constraints_.find(_label);
    if (it == constraints_.end()) {
      it = constraints_.insert(make_pair(_label, ConstraintFilter()));
    }
    it->second.addExclude(symbol);
  }

private:
  ConstraintMap constraints_;
};

typedef ConstraintFilter<fst::LogArc> LogConstraintFilter;

// ArcFilterMapper class specialization that assigns filtered weights to Zero.
template <class A, class ArcFilter>
class ArcFilterMapper {
 public:
  typedef A FromArc;
  typedef A ToArc;
  typedef ArcFilter Filter;
  typedef typename FromArc::Weight FromWeight;
  typedef typename ToArc::Weight ToWeight;

  ArcFilterMapper(const Filter &filter): filter_(filter) {};

  ToArc operator()(const FromArc &arc) const {
    if (filter_(arc)) return arc;
    else return ToArc(arc.ilabel, arc.olabel,
        ToWeight::Zero(), arc.nextstate);
  }

  fst::MapFinalAction FinalAction() const { return fst::MAP_NO_SUPERFINAL; }

  fst::MapSymbolsAction InputSymbolsAction() const { return fst::MAP_COPY_SYMBOLS; }

  fst::MapSymbolsAction OutputSymbolsAction() const { return fst::MAP_COPY_SYMBOLS;}

  uint64 Properties(uint64 props) const { return props; }

 private:
  const ArcFilter &filter_;
};

typedef ArcFilterMapper<fst::LogArc, LogConstraintFilter> ConstraintMapper;

}


#endif /* openlat_CONSTRAIN_H_ */
