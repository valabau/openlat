/*
 *   Copyright 2013, valabau
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
 * edit-distance.h
 *
 *  Created on: 03/05/2013
 *      Author: valabau
 */

#ifndef openlat_EDIT_DISTANCE_H_
#define openlat_EDIT_DISTANCE_H_

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <queue>
#include <fst/fst.h>

namespace openlat {

  template <typename T>
  size_t common_length(const vector<T>& s1, const vector<T>& s2) {
    size_t len = 0;
    const size_t lim = min(s1.size(), s2.size());
    while (len < lim) {
      if (s1[len] != s2[len]) return len;
      len++;
    }
    return len;
  }

  // structures and types for edit distance computation
  extern const char* edit_str[];
  typedef enum { EDIT_UNK = 0, EDIT_NONE, EDIT_SUB, EDIT_DEL, EDIT_INS, EDIT_NONE_NOISE, EDIT_DEL_NOISE, EDIT_INS_NOISE, MAX_EDIT} edit_t;
  extern const float default_costs[];
  extern const float default_matching_costs[];

  const float _ce = 0; // cost epsilon to penalize insertions and deletions
  //const char* word_id_str[] = { "<unk>", "<noise>", "<deleted>", "<inserted>", "<none>" };
  const char* edit_str[]        = {"UNK", "NONE", "SUB", "DEL", "INS", "NONE_NOISE", "DEL_NOISE", "INS_NOISE"};
  const float default_costs[]   = {  4,     0,      1, 1 + _ce, 1 + _ce,    0,            0,           0};
  const float default_matching_costs[] = {  4, -1,  0,     0,     0,        0,            0,           0};

  template <typename Arc>
  class EditHyp {
  public:
    Arc const * arc;
    typename Arc::StateId state;
    float cost;
    typename Arc::Weight fwd;
    typename Arc::Weight sentence_score; // score + backward_viterbi
    edit_t edit;
    std::vector<size_t> edit_count;
    size_t contiguous;
    size_t mem() { return sizeof(EditHyp<Arc>) + edit_count.size()*sizeof(size_t) + 30*sizeof(float); }
    EditHyp(): arc(0), state(fst::kNoStateId), cost(std::numeric_limits<float>::max()),
               fwd(Arc::Weight::Zero()), sentence_score(Arc::Weight::Zero()), edit(EDIT_UNK),
               edit_count(MAX_EDIT, 0), contiguous(0) {}
  };

  template <typename Arc>
  std::ostream& operator<<(std::ostream& out, const EditHyp<Arc>& val){
    out << val.cost << ":";
    out << val.sentence_score << "|";// else
    out << val.fwd << "|";
    out << std::setw(1) << edit_str[val.edit][0] << "(";
    out << std::setw(2) << ((val.arc)?val.state:fst::kNoStateId) << ")";
    out << "  ";
    return out;
  }

  template <class Arc>
  inline bool operator<(const EditHyp<Arc> &hyp1, const EditHyp<Arc> &hyp2) {
    if (hyp1.cost == hyp2.cost) {
      if (hyp1.sentence_score != hyp2.sentence_score) return hyp1.sentence_score.Value() < hyp2.sentence_score.Value();
      else return (hyp1.edit < hyp2.edit);
    }
    else return hyp1.cost < hyp2.cost;
  }

  template <typename F>
  class EditDistance {
    typedef typename F::Arc Arc;
    typedef typename Arc::Weight Weight;
    typedef typename Arc::Label Label;
    typedef typename Arc::StateId StateId;

  public:

    virtual const std::vector< EditHyp<Arc> >& getHypotheses() const = 0;

    virtual void weightsChanged() {
      backward_.clear();
      ShortestDistance(fst_, &backward_, true);
    }

    virtual void getFeatures(std::vector< std::pair< std::string, std::vector<float> > >& options) const {
      options.clear();
      const std::vector< EditHyp<Arc> >& hyps = getHypotheses();
      const fst::SymbolTable &in = *fst_.InputSymbols();
      for (StateId s = 0; s < fst_.NumStates(); s++) {
        const EditHyp<Arc>& hyp = hyps[s];
        for (fst::ArcIterator<F> ait(fst_, s); !ait.Done(); ait.Next()) {
          const Arc &arc = ait.Value();
          const Weight w = Times(hyp.fwd, Times(arc.weight, backward_[arc.nextstate]));
          options.push_back(make_pair(in.Find(arc.ilabel), getHypFeatures(w, hyp)));
        }
        // if state is final add </s>
        Weight weight = fst_.Final(s);
        if (weight != Weight::Zero()) {
          const Weight w = Times(hyp.fwd, backward_[s]);
          options.push_back(make_pair("</s>", getHypFeatures(w, hyp)));
        }
      }
    }

    virtual void updateMap(std::map<Label, std::pair<Weight, const EditHyp<Arc> *> >& words,
                   const Label &label, const Weight &weight, const EditHyp<Arc> *hyp) const
    {
      typename std::map<Label, std::pair<Weight, const EditHyp<Arc> *> >::iterator it = words.find(label);
      if (it == words.end()) {
        words.insert(make_pair(label, make_pair(weight, hyp)));
      }
      else {
        if (hyp->cost < it->second.second->cost or  (hyp->cost == it->second.second->cost and weight.Value() < it->second.first.Value())) {
          it->second.first = weight;
          it->second.second = hyp;
        }
      }
    }

    virtual void getUniqueFeatures(std::vector< std::pair< std::string, std::vector<float> > >& options) const {
      std::map<Label, std::pair<Weight, const EditHyp<Arc>*> > words;

      options.clear();
      const std::vector< EditHyp<Arc> >& hyps = getHypotheses();
      const fst::SymbolTable &in = *fst_.InputSymbols();
      Label end = in.NumSymbols();
      for (StateId s = 0; s < fst_.NumStates(); s++) {
        const EditHyp<Arc>& hyp = hyps[s];
        for (fst::ArcIterator<F> ait(fst_, s); !ait.Done(); ait.Next()) {
          const Arc &arc = ait.Value();
          const Weight w = Times(hyp.fwd, Times(arc.weight, backward_[arc.nextstate]));
          updateMap(words, arc.ilabel, w, &hyp);
        }
        // if state is final add </s>
        Weight weight = fst_.Final(s);
        if (weight != Weight::Zero()) {
          const Weight w = Times(hyp.fwd, backward_[s]);
          updateMap(words, end, w, &hyp);
        }
      }

      typename std::map<Label, std::pair<Weight, const EditHyp<Arc>*> >::const_iterator it;
      for (it = words.begin(); it != words.end(); ++it) {
        const std::string label = (it->first == end) ? "</s>" : in.Find(it->first);
        options.push_back(make_pair(label, getHypFeatures(it->second.first, *it->second.second)));
      }

    }

    virtual std::vector<float> getHypFeatures(const Weight& w, const EditHyp<Arc>& hyp) const {
      std::vector<float> features;
      features.push_back(-w.Value());
      features.push_back(hyp.contiguous);
      features.push_back(hyp.cost);
      for (size_t i = 0; i < MAX_EDIT; i++) {
        features.push_back(hyp.edit_count[i]);
      }
      return features;
    }

    virtual void appendFeatureNames(std::vector<std::string> &feature_names) const {
      feature_names.push_back("contiguous");
      feature_names.push_back("edit_distance");
      for (size_t i = 0; i < MAX_EDIT; i++) {
        std::string name = std::string("num_") + std::string(edit_str[i]);
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        feature_names.push_back(name);
      }
    }


    virtual void setPrefix(const std::vector<Label>& reference) = 0;

  EditDistance(const F& fst): fst_(fst) {
      if (not (fst.Properties(fst::kTopSorted | fst::kCoAccessible, false) & (fst::kTopSorted | fst::kCoAccessible))) {
        LOG(ERROR) << "The fst is not sorted, the edit distance algorithm will not work";
      }
      if (fst_.Start() != 0) {
        LOG(ERROR) << "The start state is not zero, the edit distance algorithm will not work";
      }
      weightsChanged();
    }

  virtual ~EditDistance() {};
  public:
    std::vector<Label> reference_;
    std::vector<Weight> backward_;
    const F& fst_;
  };

  template <>
  std::vector<float> EditDistance<LogLinearFst>::getHypFeatures(const Weight& w, const EditHyp<LogLinearArc>& hyp) const {
    std::vector<float> features;
    for (size_t i = 0; i < w.Length(); ++i) {
      features.push_back(-w.Value(i));
    }
    features.push_back(-w.Value());
    features.push_back(hyp.contiguous);
    features.push_back(hyp.cost);
    for (size_t i = 0; i < MAX_EDIT; i++) {
      features.push_back(hyp.edit_count[i]);
    }
    return features;
  }

  template <>
  void EditDistance<LogLinearFst>::appendFeatureNames(std::vector<std::string> &feature_names) const {
    feature_names.push_back("combined");
    feature_names.push_back("contiguous");
    feature_names.push_back("edit_distance");
    for (size_t i = 0; i < MAX_EDIT; i++) {
      std::string name = std::string("num_") + std::string(edit_str[i]);
      std::transform(name.begin(), name.end(), name.begin(), ::tolower);
      feature_names.push_back(name);
    }
  }


  template <typename F>
  class EditDistanceTrellis: public EditDistance<F> {
    typedef EditDistance<F> E;
    typedef typename F::Arc Arc;
    typedef typename Arc::Weight Weight;
    typedef typename Arc::Label Label;
    typedef typename Arc::StateId StateId;

  public:
    /** Prints trellis
     * @param out output stream
     */
    void printTrellis(std::ostream &out) const {
      const fst::SymbolTable &in = *E::fst_.InputSymbols();
      // print trellis_
      out.setf(std::ios::left);
      out << "---------  Trellis ---------\n";
      out << "Trellis\n\n";
      out << "Reference: size:" << E::reference_.size() << "\n";
      for (int i = (int)E::reference_.size(); i >= 0; i--) {
        out << std::setw(5) <<  ((i==0)?"":in.Find(E::reference_[i-1])) << ": ";
        for (StateId j = 0; j < E::fst_.NumStates(); j++) {
          out << trellis_[i][j].cost << ":";
          out << trellis_[i][j].sentence_score << "|";// else
          out << trellis_[i][j].fwd << "|";
          out << std::setw(1) << edit_str[trellis_[i][j].edit][0] << "(";
          out << std::setw(2) << ((trellis_[i][j].arc)?trellis_[i][j].state:fst::kNoStateId) << ")";
          out << "  ";
        }
        out << "\n";
      }
      out << "-----     ";
      for (StateId j = 0; j < E::fst_.NumStates(); j++) {
        out << std::setw(12) << j;
      }
      out << "\n";
      out << "----------------------------\n";
    }


    /** Initializes the first row of the trellis.
     */
    void initializeTrellis() {
      trellis_.reserve(1);
      trellis_.resize(1);
      /* Initialize trellis row to default values*/
      trellis_[0].clear();
      trellis_[0].reserve(E::fst_.NumStates());
      trellis_[0].resize(E::fst_.NumStates());

      std::cerr << "Trellis size at row 0*: " << (1.0*trellis_.capacity()*trellis_[0].capacity()*sizeof(EditHyp<Arc>)/(1024*1024)) << "MB" << std::endl;

      /* Set value for inital state */
      {
        EditHyp<Arc> &hyp = trellis_[0][E::fst_.Start()];
        hyp.cost  = 0;
        hyp.fwd = Weight::One();
        hyp.sentence_score = E::backward_[E::fst_.Start()];
        hyp.state = E::fst_.Start();
      }

      /* try deletions at the beginning */
      for (StateId s = 0; s < E::fst_.NumStates(); s++) {
        /* compute costs for each arc */
        const EditHyp<Arc>& prev_hyp = trellis_[0][s];
        for (fst::ArcIterator<F> ait(E::fst_, s); !ait.Done(); ait.Next()) {
          const Arc arc = ait.Value();
          /* label 0 is reserved for epsilon */
          edit_t edit_del = (arc.ilabel == 0)?EDIT_DEL_NOISE:EDIT_DEL;

          EditHyp<Arc> hyp;
          hyp.fwd = Times(prev_hyp.fwd, arc.weight);
          hyp.cost = prev_hyp.cost + default_costs[edit_del];
          hyp.edit = edit_del;
          hyp.arc = &arc;
          hyp.state = s;
          hyp.sentence_score = Times(hyp.fwd, E::backward_[arc.nextstate]);
          hyp.edit_count = prev_hyp.edit_count;
          hyp.edit_count[edit_del]++;
          hyp.contiguous = (edit_del == EDIT_DEL)? 0 : prev_hyp.contiguous;

          if (hyp < trellis_[0][arc.nextstate]) {
            trellis_[0][arc.nextstate] = hyp;
          }
        }
      }
    }

    void fillTrellisRow(size_t w) {
      trellis_.reserve(w+1);
      trellis_.resize(w+1);
      /* Initialize trellis row to default values*/
      trellis_[w].clear();
      trellis_[w].reserve(E::fst_.NumStates());
      trellis_[w].resize(E::fst_.NumStates());

      size_t ehmem = trellis_[0][0].mem();
      std::cerr << "Trellis size at row " << w << ": " << (1.0*trellis_.capacity()*trellis_[w].capacity()*ehmem/(1024*1024)) << "MB" << std::endl;

      /* try all word inserted */
      {
        edit_t edit_ins = (E::reference_[w-1] == 0)?EDIT_INS_NOISE:EDIT_INS;
        const EditHyp<Arc>& prev_hyp = trellis_[w-1][0];
        EditHyp<Arc> hyp;
        hyp.fwd = prev_hyp.fwd;
        hyp.cost = prev_hyp.cost + default_costs[edit_ins];
        hyp.edit = edit_ins;
        hyp.arc = 0;
        hyp.state = 0;
        hyp.sentence_score = Times(hyp.fwd, E::backward_[0]);
        hyp.edit_count = prev_hyp.edit_count;
        hyp.edit_count[edit_ins]++;
        hyp.contiguous = (edit_ins == EDIT_INS)? 0 : prev_hyp.contiguous;

        if (hyp < trellis_[w][0]) trellis_[w][0] = hyp;
      }

      /* fill the rest of the trellis_ */
      /* try deletions at the beginning */
      for (StateId s = 0; s < E::fst_.NumStates(); s++) {
        /* compute costs for each arc */
        for (fst::ArcIterator<F> ait(E::fst_, s); !ait.Done(); ait.Next()) {
          const Arc arc = ait.Value();

           // try insertions
            edit_t edit_ins = (E::reference_[w-1] == 0) ? EDIT_INS_NOISE : EDIT_INS;
            {
              const EditHyp<Arc>& prev_hyp = trellis_[w-1][s];
              EditHyp<Arc> hyp;
              hyp.fwd = prev_hyp.fwd;
              hyp.cost = prev_hyp.cost + default_costs[edit_ins];
              hyp.edit = edit_ins;
              hyp.arc = 0;
              hyp.state = s;
              hyp.sentence_score = Times(hyp.fwd, E::backward_[s]);
              hyp.edit_count = prev_hyp.edit_count;
              hyp.edit_count[edit_ins]++;
              hyp.contiguous = (edit_ins == EDIT_INS)? 0 : prev_hyp.contiguous;

              if (hyp < trellis_[w][s]) trellis_[w][s] = hyp;
            }


            // try deletions
            edit_t edit_del = (arc.ilabel == 0) ? EDIT_DEL_NOISE : EDIT_DEL;
            {
              const EditHyp<Arc>& prev_hyp = trellis_[w][s];
              EditHyp<Arc> hyp;
              hyp.fwd = Times(prev_hyp.fwd, arc.weight);
              hyp.cost = prev_hyp.cost + default_costs[edit_del];
              hyp.edit = edit_del;
              hyp.arc = &arc;
              hyp.state = s;
              hyp.sentence_score = Times(hyp.fwd, E::backward_[arc.nextstate]);
              hyp.edit_count = prev_hyp.edit_count;
              hyp.edit_count[edit_del]++;
              hyp.contiguous = (edit_del == EDIT_DEL)? 0 : prev_hyp.contiguous;

              if (hyp < trellis_[w][arc.nextstate]) trellis_[w][arc.nextstate] = hyp;
            }


            // try substitution
            edit_t edit_sub = (arc.ilabel == E::reference_[w-1]) ? EDIT_NONE : EDIT_SUB;
            if (edit_del == EDIT_DEL_NOISE and edit_ins == EDIT_INS_NOISE) edit_sub = EDIT_NONE_NOISE;
            if ((arc.ilabel != 0 and E::reference_[w-1] != 0) or (arc.ilabel == 0 and E::reference_[w-1] == 0)) {
              const EditHyp<Arc>& prev_hyp = trellis_[w-1][s];
              EditHyp<Arc> hyp;
              hyp.fwd = Times(prev_hyp.fwd, arc.weight);
              hyp.cost = prev_hyp.cost + default_costs[edit_sub];
              hyp.edit = edit_sub;
              hyp.arc = &arc;
              hyp.state = s;
              hyp.sentence_score = Times(hyp.fwd, E::backward_[arc.nextstate]);
              hyp.edit_count = prev_hyp.edit_count;
              hyp.edit_count[edit_sub]++;
              hyp.contiguous = (edit_sub == EDIT_SUB)? 0 : prev_hyp.contiguous;
              if (edit_sub == EDIT_NONE) hyp.contiguous++;

              if (hyp < trellis_[w][arc.nextstate]) trellis_[w][arc.nextstate] = hyp;
            }

        }
      }
    }

    virtual const std::vector< EditHyp<Arc> >& getHypotheses() const { return trellis_.back(); }

    /** Fills the trellis for edit distance between reference and lattice.
     * @param reference the prefix of the correct transcription
     */
    virtual void setPrefix(const std::vector<Label>& reference) {
      size_t len = common_length(reference, E::reference_);
      E::reference_ = reference;

      if (len == 0) initializeTrellis();
      for (size_t w = len + 1; w <= E::reference_.size(); w++) {
        fillTrellisRow(w);
      }
    }

    EditDistanceTrellis(const F& fst): EditDistance<F>(fst) {};
    virtual ~EditDistanceTrellis() {};
  protected:
    std::vector< std::vector< EditHyp<Arc> > > trellis_;
  };


  template <typename F>
  class EditDistanceRows: public EditDistance<F> {
    typedef EditDistance<F> E;
    typedef typename F::Arc Arc;
    typedef typename Arc::Weight Weight;
    typedef typename Arc::Label Label;
    typedef typename Arc::StateId StateId;

  public:

    /** Initializes the first row of the trellis.
     */
    void initializeTrellis() {
      trellis_[1].clear();
      trellis_[1].resize(E::fst_.NumStates());

      /* Set value for inital state */
      {
        EditHyp<Arc> &hyp = trellis_[1][E::fst_.Start()];
        hyp.cost  = 0;
        hyp.fwd = Weight::One();
        hyp.sentence_score = E::backward_[E::fst_.Start()];
        hyp.state = E::fst_.Start();
      }

      /* try deletions at the beginning */
      for (StateId s = 0; s < E::fst_.NumStates(); s++) {
        /* compute costs for each arc */
        const EditHyp<Arc>& prev_hyp = trellis_[1][s];
        for (fst::ArcIterator<F> ait(E::fst_, s); !ait.Done(); ait.Next()) {
          const Arc arc = ait.Value();
          /* label 0 is reserved for epsilon */
          edit_t edit_del = (arc.ilabel == 0)?EDIT_DEL_NOISE:EDIT_DEL;

          EditHyp<Arc> hyp;
          hyp.fwd = Times(prev_hyp.fwd, arc.weight);
          hyp.cost = prev_hyp.cost + default_costs[edit_del];
          hyp.edit = edit_del;
          hyp.arc = &arc;
          hyp.state = s;
          hyp.sentence_score = Times(hyp.fwd, E::backward_[arc.nextstate]);
          hyp.edit_count = prev_hyp.edit_count;
          hyp.edit_count[edit_del]++;
          hyp.contiguous = (edit_del == EDIT_DEL)? 0 : prev_hyp.contiguous;

          if (hyp < trellis_[1][arc.nextstate]) {
            trellis_[1][arc.nextstate] = hyp;
          }
        }
      }
    }

    void fillTrellisRow(size_t w) {
      trellis_[1].clear();
      trellis_[1].resize(E::fst_.NumStates());

      /* try all word inserted */
      {
        edit_t edit_ins = (E::reference_[w-1] == 0)?EDIT_INS_NOISE:EDIT_INS;
        const EditHyp<Arc>& prev_hyp = trellis_[0][0];
        EditHyp<Arc> hyp;
        hyp.fwd = prev_hyp.fwd;
        hyp.cost = prev_hyp.cost + default_costs[edit_ins];
        hyp.edit = edit_ins;
        hyp.arc = 0;
        hyp.state = 0;
        hyp.sentence_score = Times(hyp.fwd, E::backward_[0]);
        hyp.edit_count = prev_hyp.edit_count;
        hyp.edit_count[edit_ins]++;
        hyp.contiguous = (edit_ins == EDIT_INS)? 0 : prev_hyp.contiguous;

        if (hyp < trellis_[1][0]) trellis_[1][0] = hyp;
      }

      /* fill the rest of the trellis_ */
      /* try deletions at the beginning */
      for (StateId s = 0; s < E::fst_.NumStates(); s++) {
        /* compute costs for each arc */
        for (fst::ArcIterator<F> ait(E::fst_, s); !ait.Done(); ait.Next()) {
          const Arc arc = ait.Value();

           // try insertions
            edit_t edit_ins = (E::reference_[w-1] == 0) ? EDIT_INS_NOISE : EDIT_INS;
            {
              const EditHyp<Arc>& prev_hyp = trellis_[0][s];
              EditHyp<Arc> hyp;
              hyp.fwd = prev_hyp.fwd;
              hyp.cost = prev_hyp.cost + default_costs[edit_ins];
              hyp.edit = edit_ins;
              hyp.arc = 0;
              hyp.state = s;
              hyp.sentence_score = Times(hyp.fwd, E::backward_[s]);
              hyp.edit_count = prev_hyp.edit_count;
              hyp.edit_count[edit_ins]++;
              hyp.contiguous = (edit_ins == EDIT_INS)? 0 : prev_hyp.contiguous;

              if (hyp < trellis_[1][s]) trellis_[1][s] = hyp;
            }


            // try deletions
            edit_t edit_del = (arc.ilabel == 0) ? EDIT_DEL_NOISE : EDIT_DEL;
            {
              const EditHyp<Arc>& prev_hyp = trellis_[1][s];
              EditHyp<Arc> hyp;
              hyp.fwd = Times(prev_hyp.fwd, arc.weight);
              hyp.cost = prev_hyp.cost + default_costs[edit_del];
              hyp.edit = edit_del;
              hyp.arc = &arc;
              hyp.state = s;
              hyp.sentence_score = Times(hyp.fwd, E::backward_[arc.nextstate]);
              hyp.edit_count = prev_hyp.edit_count;
              hyp.edit_count[edit_del]++;
              hyp.contiguous = (edit_del == EDIT_DEL)? 0 : prev_hyp.contiguous;

              if (hyp < trellis_[1][arc.nextstate]) trellis_[1][arc.nextstate] = hyp;
            }


            // try substitution
            edit_t edit_sub = (arc.ilabel == E::reference_[w-1]) ? EDIT_NONE : EDIT_SUB;
            if (edit_del == EDIT_DEL_NOISE and edit_ins == EDIT_INS_NOISE) edit_sub = EDIT_NONE_NOISE;
            if ((arc.ilabel != 0 and E::reference_[w-1] != 0) or (arc.ilabel == 0 and E::reference_[w-1] == 0)) {
              const EditHyp<Arc>& prev_hyp = trellis_[0][s];
              EditHyp<Arc> hyp;
              hyp.fwd = Times(prev_hyp.fwd, arc.weight);
              hyp.cost = prev_hyp.cost + default_costs[edit_sub];
              hyp.edit = edit_sub;
              hyp.arc = &arc;
              hyp.state = s;
              hyp.sentence_score = Times(hyp.fwd, E::backward_[arc.nextstate]);
              hyp.edit_count = prev_hyp.edit_count;
              hyp.edit_count[edit_sub]++;
              hyp.contiguous = (edit_sub == EDIT_SUB)? 0 : prev_hyp.contiguous;
              if (edit_sub == EDIT_NONE) hyp.contiguous++;

              if (hyp < trellis_[1][arc.nextstate]) trellis_[1][arc.nextstate] = hyp;
            }

        }
      }
    }

    virtual const std::vector< EditHyp<Arc> >& getHypotheses() const { return trellis_.back(); }

    /** Fills the trellis for edit distance between reference and lattice.
     * @param reference the prefix of the correct transcription
     */
    virtual void setPrefix(const std::vector<Label>& reference) {
      size_t len = common_length(reference, E::reference_);

      size_t w = 1;

      std::cerr << "len = " << len << " ref = " << reference.size() << "\n";

      if (len == E::reference_.size() - 1) {
        std::cerr << "len == E::reference_.size() - 1\n";
        w = len + 1;
        std::swap(trellis_[0], trellis_[1]);
      }
      else if (len == E::reference_.size() and not reference.empty()) {
        std::cerr << "len == E::reference_.size() and not E::reference_.empty()\n";
        w = len + 1;
      }
      else {
        std::cerr << "else initialize\n";
        initializeTrellis();
        w = 1;
      }

      std::cerr << "t[0][0] " << trellis_[0][0] << "\n";
      std::cerr << "t[1][0] " << trellis_[1][0] << "\n";

      E::reference_ = reference;
      for (; w <= E::reference_.size(); w++) {
        std::cerr << "w = " << w << std::endl;
        swap(trellis_[0], trellis_[1]);
        fillTrellisRow(w);
      }

    }

    EditDistanceRows(const F& fst): EditDistance<F>(fst), n_rows_(2) {
      trellis_.reserve(n_rows_);
      trellis_.resize (n_rows_);

      for (size_t i = 0; i < n_rows_; i++) {
        trellis_[i].reserve(E::fst_.NumStates());
        trellis_[i].resize(E::fst_.NumStates());
      }

    };
    virtual ~EditDistanceRows() {};
  protected:
    std::vector< std::vector< EditHyp<Arc> > > trellis_;
    size_t n_rows_;
  };

}


#endif /* openlat_EDIT_DISTANCE_H_ */
