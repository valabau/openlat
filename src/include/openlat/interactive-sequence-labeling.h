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
 * interactive-sequence-labeling.h
 *
 *  Created on: 06/03/2012
 *      Author: valabau
 */

#ifndef openlat_INTERACTIVE_SEQUENCE_LABELING_H_
#define openlat_INTERACTIVE_SEQUENCE_LABELING_H_

#include <tr1/unordered_map>
#include <set>
#include <fst/shortest-path.h>
#include <openlat/compat.h>
#include <openlat/utils.h>
#include <openlat/rmarc.h>
#include <openlat/verify.h>
#include <openlat/query.h>

namespace openlat {

template<class Arc>
void GetOutputSymbols(const fst::Fst<Arc> &fst, std::vector<std::set<typename Arc::Label> > *symbols) {
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId StateId;

  if (not fst.Properties(fst::kTopSorted, true)) {
    LOG(FATAL) << "VerifySequenceLabeling: the fst must be topologically sorted to verify sequence labeling\n";
  }

  if (not fst.Properties(fst::kAcyclic, true)) {
    LOG(FATAL) << "VerifySequenceLabeling: the fst must be acyclic to verify sequence labeling\n";
  }

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
       siter.Next())
  {
    StateId s = siter.Value();
    size_t this_length = state_length[s];

    if (this_length != static_cast<size_t>(-1)) {
      symbols->resize(this_length + 1);
      size_t na = 0;

      for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s);
           !aiter.Done();
           aiter.Next())
      {
        const Arc &arc = aiter.Value();

        (*symbols)[this_length].insert(arc.olabel);

        if (state_length[arc.nextstate] == static_cast<size_t>(-1)) state_length[arc.nextstate] = this_length + 1;
        else if (state_length[arc.nextstate] != this_length + 1) {
          LOG(FATAL) << "VerifySequenceLabeling: Fst state " << state_length[arc.nextstate]
                     << " length from " << s << " is different from a previous state";
        }

        ++na;
      }
      if (!fst.Final(s).Member()) {
        LOG(FATAL) << "VerifySequenceLabeling: Fst final weight of state " << s << " is invalid";
      }
      if (fst.Final(s).Value() < Weight::Zero().Value()) {
        if (final_length == static_cast<size_t>(-1)) final_length = this_length;
        else if (this_length != final_length) {
          LOG(FATAL) << "VerifySequenceLabeling: Fst state " << s << " length is different from a previous final state length";
        }
      }
    }
  }
}


typedef struct _sPoolV {
  VLabel hyp;
  float score;
  _sPoolV(VLabel _hyp, float _score): hyp(_hyp), score(_score) {}
  bool operator<(const _sPoolV &other) const { return score < other.score; }
} PoolV;

typedef std::tr1::unordered_map<SampleLabel, PoolV> Pool;

}

namespace std {
  template<> struct less<openlat::Pool::iterator> {
    size_t operator()(const openlat::Pool::iterator& __x, const openlat::Pool::iterator& __y) const {
      return __x->second < __y->second;
    }
  };
}

namespace openlat {

/**
 * @param clock1 first clock
 * @param clock2 second clock
 * @return return the lapse of time in ms between clock1 and clock2
 */
template <typename T>
double diffclock(T clock1, T clock2) {
  double diffticks = static_cast<double>(clock1) - static_cast<double>(clock2);
  double diffms = (diffticks * 1000) / CLOCKS_PER_SEC;
  return diffms;
}


/** @class ISystem
 *
 * @brief An interactive sequence labeling system to be queried by the oracle.
 *
 */
class ISystem {
public:
  /** Actively select a query (label, hypothesis) to be validated by the Oracle */
  virtual Query askLabel() = 0;
  /** Assign a hypothesis to a label (query), which was erroneously recognized, and return the most likely hypothesis with that constraint */
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) = 0;
  /** Recompute score for given query and return the most likely hyp and the scores of the elements */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) = 0;
  /** Accept quiven query as correct */
  virtual void acceptLabel(const Query &query) = 0;
  /** Update a hypothesis to a label (query) and return the most likely hypothesis with that constraint */
  virtual void update(const Query &query, vector<VLabel> &hyp) = 0;
  virtual ~ISystem() {};
};


/** @class Oracle
 *
 *  @brief the Oracle evaluates interactive systems that implement the \ref ISystem interface.
 */
class Oracle {
public:
  const vector<vector<VLabel> >& refs; /**< correct answers for each sample */
  vector<size_t> n_errors; /**< number of errors for each sample */
  size_t c; /**< number of corrections made */
  size_t e_c; /**< number of label errors */
  size_t e_k; /**< number of @e sentence errors */
  clock_t start_time; /**< start time of the evaluation */
  clock_t end_time; /**< end time of the evaluation */

  Oracle(const vector<vector<VLabel> >& _refs):
          refs(_refs), n_errors(refs.size()), c(0)
  {
    for (size_t i = 0; i < refs.size(); i++) n_errors[i] = refs[i].size();
    recomputeErrors();
  }

  /**< recompute errors from n_errors */
  void recomputeErrors() {
    e_c = 0;
    e_k = 0;
    for (size_t i = 0; i < n_errors.size(); i++) {
      if (n_errors[i] > 0) e_k++;
      e_c += n_errors[i];
    }
  }

  /** print error rates and other statistics */
  void printStats(size_t s) {
    recomputeErrors();
    end_time = clock();
    cout << s << " " << c << " " << e_c << " " << e_k << " " << diffclock(end_time, start_time)/1000 << endl;
  }

  /** evaluate an interactive system */
  void evaluate(ISystem *system) {
    start_time = clock();

    // compute initial error rates and print initial statistics
    for (size_t s = 0; s < refs.size(); s++) {
      Query query(s, fst::SymbolTable::kNoSymbol, fst::SymbolTable::kNoSymbol);
      vector<VLabel> hyp;
      system->fixLabel(query, hyp);
      n_errors[query.label.sample] = edit_distance(refs[query.label.sample], hyp);
    }
    printStats(0);

    size_t total_labels = 0; // total labels to be queried
    for (size_t i = 0; i < refs.size(); i++) total_labels += refs[i].size();

    // start interactive labeling until total labels are queried
    for (size_t s = 0; s < total_labels; s++) {
      // ask the system for the next query
      Query query = system->askLabel();
      // cerr << "query: " << query.label.sample << " " << query.label.label << " " << query.hyp << " is " << refs[query.label.sample][query.label.label] << "\n";

      // if the query is wrong
      if (query.hyp != refs[query.label.sample][query.label.label]) {
        c++; // increment the number of corrections made
        query.hyp = refs[query.label.sample][query.label.label]; // assign the correct hyp
        vector<VLabel> hyp;
        system->fixLabel(query, hyp); // fix the hypothesis in the system
        // recompute the number of errors for the given sample
        n_errors[query.label.sample] = edit_distance(refs[query.label.sample], hyp);
      }
      // if the query is right then accept the query
      else {
//        system->acceptLabel(query);
        vector<VLabel> hyp;
        system->fixLabel(query, hyp); // this is necessary to acknowledge the system
      }
      // print the updated stats
      printStats(s + 1);
    }
  }
};



/** @class System
 *  @brief Base class for interactive systems
 */
template <class Arc, class Filter>
class System: public ISystem {
public:
  vector<fst::VectorFst<Arc> *> fsts; //< store all fsts samples
  unordered_map<SampleLabel, PoolV> elems; //< pool of hypotheses
  vector< vector<bool> > assigned_labels; //< store already assigned labels for each sample

  System(vector<fst::VectorFst<Arc> *> _fsts): fsts(_fsts) {
    assigned_labels.resize(fsts.size());
    for (size_t s = 0; s < fsts.size(); s++) {
      // obtain the number of labels for this sample
      size_t n_labels = VerifySequenceLabeling(*fsts[s]);
      // allocate space for the output labels and the pool of hypotheses
      assigned_labels[s].resize(n_labels);
      for (size_t l = 0; l < n_labels; l++) {
        elems.insert(make_pair(SampleLabel(s, l), PoolV(0, 0)));
      }
    }

  }

  virtual Query askLabel() = 0;
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) = 0;
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) = 0;

  /** accept query as a correct solution */
  virtual void acceptLabel(const Query &query) {
    // mark current label as validated
    assigned_labels[query.label.sample][query.label.label] = true;

    // retrieve the fst and remove all hypotheses not compatible with the current query
    fst::VectorFst<Arc> * fst = fsts[query.label.sample];
    Filter filter(query.label.label, query.hyp);
    RmArc<Arc, Filter>(fst, RmArcOptions<Arc, Filter>(filter));
    assert_bt(fst->NumStates() > 0, "ERROR: empty lattice after accepting label\n");
  }

  /** update the system with the query and return the updated hypothesis */
  virtual void update(const Query &query, vector<VLabel> &hyp) {
    // accept if the query is known
    if (query.label.label != fst::kNoLabel) {
      acceptLabel(query);
    }

    // recompute system scores and update the pool of hypotheses
    vector<float> scores;
    recomputeScores(query, hyp, scores);
    for (size_t i = 0; i < hyp.size(); i++) {
      SampleLabel label(query.label.sample, i);
      Pool::iterator it = elems.find(label);
      it->second.hyp = hyp[i];
      it->second.score = scores[i];
    }
  }
  virtual ~System() {};
};


/** @class LocalSystem
 *  @brief A local system selects the query based exclusively on the current sample
 *
 *  The local system iterates over the samples making just one query per sample.
 *  In the next iteration it performs another query per sample, and so on.
 */
template <class Arc, class Filter, typename Recompute>
class LocalSystem: public System<Arc, Filter> {
  typedef System<Arc, Filter> Base;
public:
  LocalSystem(vector<fst::VectorFst<Arc> *>& _fsts): Base(_fsts), next_sample(0), pool(Base::fsts.size()) {
    for (Pool::iterator it = Base::elems.begin(); it != Base::elems.end(); ++it) {
      pool[it->first.sample].push_back(it);
    }
  }

  /** Ask next label */
  virtual Query askLabel() {
    Query query;

    // process next sample in the queue
    query.label.sample = next_sample;
    next_sample = (next_sample + 1) % Base::fsts.size();

    // pop min element from pool[query.sample]
    pop_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), std::less<Pool::iterator>());

    // return query and remove it from the queue
    query.label.label = pool[query.label.sample].back()->first.label;
    query.hyp         = pool[query.label.sample].back()->second.hyp;
    pool[query.label.sample].pop_back();
    return query;
  }

  /** Assign query and update re-sort the queue jut for the current sample */
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    Base::update(query, hyp);
    make_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), std::less<Pool::iterator>());
  }

  /** Recompute scores using the Recompute template argument */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    fst::VectorFst<Arc> *fst = Base::fsts[query.label.sample];
    // create a dummy dead state for the mappers to disable an arc
    typename fst::VectorFst<Arc>::StateId dead_state = fst->AddState();
    Recompute()(*fst, dead_state,
        Base::assigned_labels[query.label.sample], hyp, scores);
    // delete the dummy state
    fst->DeleteStates(vector<typename fst::VectorFst<Arc>::StateId>(1, dead_state));
  }

  virtual ~LocalSystem() {};

protected:
  size_t next_sample; //< which is the next sample to process
  vector< vector<Pool::iterator>  > pool; //< a pool of queries for each sample
};

/** @class LocalSystem
 *  @brief A global system selects the query based on the whole set of samples
 *
 */
template <class Arc, class Filter, typename Recompute>
class GlobalSystem: public System<Arc, Filter> {
public:
  typedef System<Arc, Filter> Base;
  GlobalSystem(std::vector<fst::VectorFst<Arc> *>& _fsts): Base(_fsts) {
    for (Pool::iterator it = Base::elems.begin(); it != Base::elems.end(); ++it) {
      pool.push_back(it);
    }
  }

  /** Ask next label */
  virtual Query askLabel() {
    Query query;
    // pop min element from the global pool
    pop_heap(pool.begin(), pool.end(), std::less<Pool::iterator>());
    // return query
    query.label.sample = pool.back()->first.sample;
    query.label.label  = pool.back()->first.label;
    query.hyp          = pool.back()->second.hyp;
    pool.pop_back();
    return query;
  }


  /** Assign query and update re-sort the queue jut for the current sample */
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    Base::update(query, hyp);
    make_heap(pool.begin(), pool.end(), std::less<Pool::iterator>());
  }

  /** Recompute scores using the Recompute template argument */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    fst::VectorFst<Arc> *fst = Base::fsts[query.label.sample];
    // create a dummy dead state for the mappers to disable an arc
    typename fst::VectorFst<Arc>::StateId dead_state = fst->AddState();
    Recompute()(*fst, dead_state,
        Base::assigned_labels[query.label.sample], hyp, scores);
    // delete the dummy state
    fst->DeleteStates(vector<typename fst::VectorFst<Arc>::StateId>(1, dead_state));
  }
  virtual ~GlobalSystem() {};

protected:
  vector<Pool::iterator> pool; //< a global pool of queries
};

template <class Arc, class Filter>
struct RecomputeRandom {
  void operator()(const fst::VectorFst<Arc> &fst, typename fst::VectorFst<Arc>::StateId, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    ShortestPath(fst, hyp, scores);
    for (size_t l = 0; l < scores.size(); l++) {
      scores[l] = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
    }
  }
};

template <class Arc, class Filter>
struct RecomputeGreedy {
  void operator()(const fst::VectorFst<Arc> &fst, typename fst::VectorFst<Arc>::StateId, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    ShortestPath(fst, hyp, scores);
  }
};

template <class Arc, class Filter>
struct RecomputeSequential {
  void operator()(const fst::VectorFst<Arc> &fst, typename fst::VectorFst<Arc>::StateId, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    ShortestPath(fst, hyp, scores);
    for (size_t l = 0; l < hyp.size(); l++) {
      scores[l] = -l;
    }
  }
};

template <class Arc, class Filter>
struct RecomputeExpectation {
  void operator()(const fst::VectorFst<Arc> &fst, typename fst::VectorFst<Arc>::StateId dead_state, const vector<bool>& assigned_labels, vector<VLabel> &hyp, vector<float> &scores) {
    typedef FilterMapper<Arc, Filter> Mapper;

    scores.clear();
    scores.resize(assigned_labels.size(), Arc::Weight::Zero().Value());

    hyp.resize(assigned_labels.size());

    std::vector< std::set<typename Arc::Label> > members;
    GetOutputSymbols(fst, &members);

  //  cout << "assigned labels " << assigned_labels.to_string() << "\n";

    float max_score = Arc::Weight::Zero().Value();
    for (size_t label = 0; label < assigned_labels.size(); label++) {
      vector<VLabel> best_hyp(assigned_labels.size(), fst::SymbolTable::kNoSymbol);

      if (not assigned_labels[label]) {
        float max_img_score = Arc::Weight::Zero().Value();
        for (typename std::set<typename Arc::Label>::const_iterator mem_it = members[label].begin();
            mem_it != members[label].end(); ++mem_it)
        {
          size_t member = *mem_it;
          Filter filter(label, member);
          Mapper mapper(filter, dead_state);
          fst::MapFst<Arc, Arc, Mapper> _fst(fst, mapper);

          vector<VLabel> lhyp;
          vector<float> lscores;
          float score = ShortestPath(_fst, lhyp, lscores);
  //            if (label == 0 and score < 1e-5) {
  //              fst_kbest->write(cout);
  //            }
  //            cout << "count " << count << "\n";
  //            cout << "score l" << label << " i" << member << ": " << score << " ";
  //            cout << sentence_to_string(lhyp);
  //            cout << "\n";

  //            cout << "score l" << label << " i" << member << " = " << score << "\n";
          scores[label] = add_log(scores[label], score);
          if (score < max_img_score) {
            max_img_score = score;
            copy(lhyp.begin(), lhyp.end(), best_hyp.begin());
          }
        }

  //      cout << "max score l" << label << ": " << scores[label] << " ";
  //      cout << sentence_to_string(best_hyp);
  //      cout << "\n";
        if (scores[label] < max_score) {
          max_score = scores[label];
          copy(best_hyp.begin(), best_hyp.end(), hyp.begin());
        }
      } // if not assigned label
    }

  //  cout << "max: " << max_score << "\n";
  //  cout << "scores: ";
  //  for (size_t s = 0; s < scores.size(); s++) {
  //    cout << " " << scores[s];
  //  }
  //  cout << "\n";

    vector<float> tscores;
    ShortestPath(fst, hyp, tscores);
  }
};


}

#endif /* openlat_INTERACTIVE_SEQUENCE_LABELING_H_ */
