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
#include <openlat/constrain.h>

namespace openlat {

/** Obtain a set of symbols from the fst
 * @param[in] fst a fst
 * @param[out] a set of symbols
 */
template<class Arc>
void GetOutputSymbolsOld(const fst::Fst<Arc> &fst, std::vector<std::set<typename Arc::Label> > *symbols) {
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId StateId;

  if (not fst.Properties(fst::kTopSorted, true)) {
    LOG(FATAL) << "VerifySequenceLabeling: the fst must be topologically sorted to verify sequence labeling\n";
  }

  if (not fst.Properties(fst::kAcyclic, true)) {
    LOG(FATAL) << "VerifySequenceLabeling: the fst must be acyclic to verify sequence labeling\n";
  }

  symbols->clear();

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

        if (arc.weight != Weight::Zero()) {
          (*symbols)[this_length].insert(arc.olabel);

          if (state_length[arc.nextstate] == static_cast<size_t>(-1)) state_length[arc.nextstate] = this_length + 1;
          else if (state_length[arc.nextstate] != this_length + 1) {
            LOG(FATAL) << "VerifySequenceLabeling: Fst state " << state_length[arc.nextstate]
                       << " length from " << s << " is different from a previous state";
          }
        }

        ++na;
      }
      if (!fst.Final(s).Member()) {
        LOG(FATAL) << "VerifySequenceLabeling: Fst final weight of state " << s << " is invalid";
      }
      if (fst.Final(s) != Weight::Zero()) {
        if (final_length == static_cast<size_t>(-1)) final_length = this_length;
        else if (this_length != final_length) {
          LOG(FATAL) << "VerifySequenceLabeling: Fst state " << s << " length is different from a previous final state length";
        }
      }
    }
  }
}

template<class Arc>
void GetOutputSymbols(const fst::Fst<Arc> &fst, std::vector<std::set<typename Arc::Label> > *symbols) {
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId StateId;

  symbols->clear();

  for (fst::StateIterator< fst::Fst<Arc> > siter(fst);
       !siter.Done();
       siter.Next())
  {
    for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, siter.Value());
         !aiter.Done();
         aiter.Next())
    {
      const Arc &arc = aiter.Value();

      if (arc.ilabel != fst::kNoLabel and arc.ilabel > 0) {

        if (arc.weight != Weight::Zero()) {
          size_t label = arc.ilabel - 1;
          if ((*symbols).size() <= label) (*symbols).resize(label + 1);
          (*symbols)[label].insert(arc.olabel);
        }
      }
    }
  }
}


/**
 * structure to manage the pool of hypothesis
 */
typedef struct _sPoolV {
  VLabel hyp;
  float score; //< score of the hypothesis in the log semiring (the smaller the better)
  _sPoolV(VLabel _hyp, float _score): hyp(_hyp), score(_score) {}
  bool operator<(const _sPoolV &other) const { return score < other.score; }
  bool operator>(const _sPoolV &other) const { return score > other.score; }
} PoolV;

typedef std::tr1::unordered_map<SampleLabel, PoolV> Pool;

/**
 * structure function to sort the Pool acording to label
 */
struct sort_pool_by_label {
  size_t operator()(const Pool::const_iterator& __x, const Pool::const_iterator& __y) const {
    return __x->first.label > __y->first.label;
  }
};

/**
 * structure function to sort the Pool acording to label
 */
struct sort_pool_by_label_reverse {
  size_t operator()(const Pool::const_iterator& __x, const Pool::const_iterator& __y) const {
    return __x->first.label < __y->first.label;
  }
};

/**
 * structure function to sort the Pool
 */
struct sort_pool_by_score {
  size_t operator()(const Pool::const_iterator& __x, const Pool::const_iterator& __y) const {
    return __x->second > __y->second;
  }
};

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
  /** Accept given query as correct */
  virtual void acceptLabel(const Query &query) = 0;
  /** Update a hypothesis to a label (query) and return the most likely hypothesis with that constraint */
  virtual void update(const Query &query, vector<VLabel> &hyp) = 0;
  /** Recompute score for given query and return the most likely hyp and the scores of the elements */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) = 0;
  virtual ~ISystem() {};
};


/** @class Oracle
 *
 *  @brief the Oracle evaluates interactive systems that implement the \ref ISystem interface.
 */
class Oracle {
public:
  const std::vector<std::vector<VLabel> >& refs; /**< correct answers for each sample */
  const fst::SymbolTable *isymbols; /**< input symbol table */
  const fst::SymbolTable *osymbols; /**< output symbol table */
  std::vector<size_t> n_errors; /**< number of errors for each sample */
  size_t c; /**< number of corrections made */
  size_t n_labels; /**< total number of labels */
  size_t e_c; /**< number of label errors */
  size_t e_k; /**< number of @e sentence errors */
  clock_t start_time; /**< start time of the evaluation */
  clock_t end_time;   /**< end time of the evaluation */

  Oracle(const std::vector<std::vector<VLabel> >& _refs,
         const fst::SymbolTable *_isymbols,
         const fst::SymbolTable *_osymbols):
          refs(_refs), isymbols(_isymbols), osymbols(_osymbols),
          n_errors(refs.size()), c(0), n_labels(0),
          e_c(0), e_k(0), start_time(0), end_time(0)
  {
    // initially the error is 100%
    for (size_t i = 0; i < refs.size(); i++) {
      n_labels += refs[i].size();
      n_errors[i] = refs[i].size();
    }
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

  float error(int count, int total) { return 100.0*static_cast<float>(count)/static_cast<float>(total); }

  /** print error rates and other statistics */
  void printStats(size_t s, const Query& query) {
    int pos = (s == 0)?-1:((s-1) / refs.size());
    int ls = (query.label.sample == size_t(-1))?-1:query.label.sample;

    recomputeErrors();
    end_time = clock();
    std::cerr << pos;
    std::cerr << " " << error(s, n_labels)      << " (" << s   << ")";
    std::cerr << " Q(s_" << ls << "[";
    if (isymbols == 0)  std::cerr << query.label.label ;
    else std::cerr << isymbols->Find(query.label.label+1);
    std::cerr << "]==";
    if (osymbols == 0)  std::cerr << query.hyp ;
    else std::cerr << osymbols->Find(query.hyp);
    std::cerr <<  ")";
    std::cerr << " " << error(c, n_labels)      << " (" << c   << ")";
    std::cerr << " " << error(e_c, n_labels)    << " (" << e_c << ")";
    std::cerr << " " << error(e_k, refs.size()) << " (" << e_k << ")";
    std::cerr << " " << diffclock(end_time, start_time)/1000 << std::endl;

    std::map<size_t, size_t> histogram;
    for (size_t i = 0; i < n_errors.size(); i++) {
      histogram[n_errors[i]]++;
    }
    std::cerr << "{";
    for (std::map<size_t, size_t>::const_iterator it = histogram.begin(); it != histogram.end(); ++it) {
      std::cerr << " " << it->first << ":" << it->second << ",";
    }
    std::cerr << "}" << std::endl;

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
//      if (n_errors[query.label.sample] > 0) cout << "filt " << query.label.sample << " " << n_errors[query.label.sample] << "\n";
    }
    printStats(0, Query());

    // start interactive labeling until total labels are queried
    for (size_t s = 0; s < n_labels; s++) {
      // ask the system for the next query
      Query query = system->askLabel();
      // cerr << "query: " << query.label.sample << " " << query.label.label << " " << query.hyp << " is " << refs[query.label.sample][query.label.label] << "\n";

      // if the query is wrong
      if (query.hyp != refs[query.label.sample][query.label.label]) {
//        const vector<VLabel> &ref = refs[query.label.sample];
        c++; // increment the number of corrections made
        query.hyp = refs[query.label.sample][query.label.label]; // assign the correct hyp
        vector<VLabel> hyp;
        system->fixLabel(query, hyp); // fix the hypothesis in the system
        // recompute the number of errors for the given sample
        size_t dist = edit_distance(refs[query.label.sample], hyp);
//        size_t n_err = n_errors[query.label.sample];
        n_errors[query.label.sample] = dist;
      }
      // if the query is right then accept the query
      else {
//        system->acceptLabel(query);
        vector<VLabel> hyp;
        system->fixLabel(query, hyp); // this is necessary to acknowledge the system
      }
      // print the updated stats
      printStats(s + 1, query);
    }
  }
};



/** @class System
 *  @brief Base class for interactive systems
 */
template <class Arc, class Filter, typename Recompute>
class System: public ISystem {
public:
  std::vector<fst::MutableFst<Arc> *> fsts; //< store all fsts samples
  std::vector<Filter> constraints; //< store all fsts constraints
  std::tr1::unordered_map<SampleLabel, PoolV> elems; //< pool of hypotheses
  std::vector< std::vector<bool> > assigned_labels; //< store already assigned labels for each sample
  Recompute recompute;

  System(const std::vector<const fst::MutableFst<Arc> *> &_fsts) {
    for (size_t i = 0; i < _fsts.size(); i++) {
      fsts.push_back(recompute.prepare(*_fsts[i]));
    }
    constraints.resize(fsts.size());
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

  virtual ~System() {
    for (size_t i = 0; i < fsts.size(); i++) {
      delete fsts[i];
    }
  };

  virtual Query askLabel() = 0;
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) = 0;

  /** Recompute scores using the Recompute template argument */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    recompute(*fsts[query.label.sample], constraints[query.label.sample], assigned_labels[query.label.sample], hyp, scores);
  }

  /** accept query as a correct solution */
  virtual void acceptLabel(const Query &query) {
    // mark current label as validated
    assigned_labels[query.label.sample][query.label.label] = true;

    // retrieve the fst and remove all hypotheses not compatible with the current query
    constraints[query.label.sample].setMatch(query.label.label, query.hyp);
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
};


/** @class LocalSystem
 *  @brief A local system selects the query based exclusively on the current sample
 *
 *  The local system iterates over the samples making just one query per sample.
 *  In the next iteration it performs another query per sample, and so on.
 */
template <class Arc, class Filter, typename Recompute, typename Less = sort_pool_by_score>
class LocalSystem: public System<Arc, Filter, Recompute> {
  typedef System<Arc, Filter, Recompute> Base;
public:
  LocalSystem(std::vector<const fst::MutableFst<Arc> *>& _fsts): Base(_fsts), next_sample(0), pool(Base::fsts.size()) {
    for (Pool::const_iterator it = Base::elems.begin(); it != Base::elems.end(); ++it) {
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
    pop_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), Less());

    // return query and remove it from the queue
    query.label.label = (pool[query.label.sample].back())->first.label;
    query.hyp         = (pool[query.label.sample].back())->second.hyp;
    pool[query.label.sample].pop_back();
    return query;
  }

  /** Assign query and update re-sort the queue jut for the current sample */
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    Base::update(query, hyp);
    std::vector<Pool::const_iterator> &ps = pool[query.label.sample];
    make_heap(ps.begin(), ps.end(), Less());
  }

  virtual ~LocalSystem() {};

protected:
  size_t next_sample; //< which is the next sample to process
  vector< vector<Pool::const_iterator>  > pool; //< a pool of queries for each sample
};

/** @class GlobalSystem
 *  @brief A global system selects the query based on the whole set of samples
 *
 */
template <class Arc, class Filter, typename Recompute, typename Less = sort_pool_by_score>
class GlobalSystem: public System<Arc, Filter, Recompute> {
public:
  typedef System<Arc, Filter, Recompute> Base;
  GlobalSystem(std::vector<const fst::MutableFst<Arc> *>& _fsts): Base(_fsts) {
    for (Pool::const_iterator it = Base::elems.begin(); it != Base::elems.end(); ++it) {
      pool.push_back(it);
    }
  }

  /** Ask next label */
  virtual Query askLabel() {
    Query query;
    // pop min element from the global pool
    pop_heap(pool.begin(), pool.end(), Less());
    // return query
    query.label.sample = (pool.back())->first.sample;
    query.label.label  = (pool.back())->first.label;
    query.hyp          = (pool.back())->second.hyp;
    pool.pop_back();
    return query;
  }


  /** Assign query and update re-sort the queue jut for the current sample */
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    Base::update(query, hyp);
    make_heap(pool.begin(), pool.end(), Less());
  }
  virtual ~GlobalSystem() {};

protected:
  vector<Pool::const_iterator> pool; //< a global pool of queries
};

template <class Arc, class Filter>
struct RecomputeRandom {
  fst::MutableFst<Arc> *prepare(const fst::MutableFst<Arc> &fst) { return fst.Copy(); }
  void operator()(const fst::Fst<Arc> &fst, const Filter &filter, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    typedef ArcFilterMapper<Arc, Filter> Mapper;
    Mapper mapper(filter);
    fst::MapFst<Arc, Arc, Mapper> _fst(fst, mapper);

    ShortestPath(_fst, hyp, scores);
    for (size_t l = 0; l < scores.size(); l++) {
      scores[l] = -static_cast<float>(rand())/static_cast<float>(RAND_MAX);
    }
  }
};

template <class Arc, class Filter>
struct RecomputeGreedy {
  fst::MutableFst<Arc> *prepare(const fst::MutableFst<Arc> &fst) { return fst.Copy(); }
  void operator()(const fst::Fst<Arc> &fst, const Filter &filter, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    typedef ArcFilterMapper<Arc, Filter> Mapper;
    Mapper mapper(filter);
    fst::MapFst<Arc, Arc, Mapper> _fst(fst, mapper);

    ShortestPath(_fst, hyp, scores);
  }
};

template <class Arc, class Filter>
struct RecomputeSequential {
  fst::MutableFst<Arc> *prepare(const fst::MutableFst<Arc> &fst) { return fst.Copy(); }
  void operator()(const fst::Fst<Arc> &fst, const Filter &filter, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    typedef ArcFilterMapper<Arc, Filter> Mapper;
    Mapper mapper(filter);
    fst::MapFst<Arc, Arc, Mapper> _fst(fst, mapper);

    ShortestPath(_fst, hyp, scores);
    for (size_t l = 0; l < hyp.size(); l++) {
      scores[l] = -scores[l];
    }
  }
};


template <class Arc>
typename Arc::Weight FstMass(const fst::Fst<Arc> &fst) {
  std::vector<typename Arc::Weight> distance;
  fst::ShortestDistance(fst, &distance, true);
  return distance[fst.Start()];
}

template <class Arc, class Filter, bool reverse = false>
struct RecomputeSequentialExpectation {
  fst::MutableFst<Arc> *prepare(const fst::MutableFst<Arc> &fst) {
    if (reverse) {
      // first reverse the original fst
      fst::MutableFst<Arc> *rfst = new fst::VectorFst<Arc>();
      fst::Reverse(fst, rfst);

      // now normalize
      fst::MutableFst<Arc> *nfst = new fst::VectorFst<Arc>();
      DeterminizeAndNormalize(*rfst, nfst);
      delete rfst;
      return nfst;
    }
    else {
      fst::MutableFst<Arc> *nfst = new fst::VectorFst<Arc>();
      DeterminizeAndNormalize(fst, nfst);
      return nfst;
    }
  }
  void operator()(const fst::Fst<Arc> &fst, const Filter &filter, const vector<bool>& assigned_labels, vector<VLabel> &hyp, vector<float> &scores) {
    typedef ArcFilterMapper<Arc, Filter> Mapper;
    typedef typename Arc::Weight Weight;
    typedef typename Arc::Label Label;
    typedef typename Arc::StateId StateId;

    // reset scores and hyp
    scores.clear();
    scores.resize(assigned_labels.size(), -Arc::Weight::Zero().Value());
    hyp.resize(assigned_labels.size());

    StateId current_state = fst.Start();
    for (size_t l = 0; l < assigned_labels.size(); l++) {
      size_t label = l;
      if (reverse) label = assigned_labels.size() - l - 1;

      Weight max_weight = Weight::Zero();

      for (fst::ArcIterator < fst::Fst<Arc> > aiter(fst, current_state); !aiter.Done(); aiter.Next()) {
        const Arc &arc = aiter.Value();
        Label assigned = filter.getMatch(label);
        if (assigned_labels[label]) {
          assert_bt(assigned != fst::kNoLabel, "Unexpected unassigned label");
          if (arc.olabel == assigned) {
            max_weight = arc.weight;
            scores[label] = -arc.weight.Value();
            hyp[label] = arc.olabel;
            current_state = arc.nextstate;
          }
        }
        else if (arc.weight != Weight::Zero()) {
//          std::string member_str = fst.OutputSymbols()->Find(arc.olabel);
//          std::string label_str = fst.InputSymbols()->Find(arc.ilabel);

          if (arc.weight.Value() <= max_weight.Value()) {
            max_weight = arc.weight;
            scores[label] = -arc.weight.Value();
            hyp[label] = arc.olabel;
            current_state = arc.nextstate;
          }
        }
      }
    } // if not assigned label

//    std::cerr << "output: " << syms_to_string(hyp, fst.OutputSymbols()) << std::endl;
  }
};


template <class Arc, class Filter>
struct RecomputeExpectation {
  fst::MutableFst<Arc> *prepare(const fst::MutableFst<Arc> &fst) { return fst.Copy(); }
  void operator()(const fst::Fst<Arc> &fst, const Filter &filter, const vector<bool>& assigned_labels, vector<VLabel> &hyp, vector<float> &scores) {
    typedef FilterMapper<Arc, Filter> Mapper;

    // reset scores and hyp
    scores.clear();
    scores.resize(assigned_labels.size(), -Arc::Weight::Zero().Value());
    hyp.resize(assigned_labels.size());
    float max_score = -Arc::Weight::Zero().Value();

    // obtain the available members for each label
    std::vector< std::set<typename Arc::Label> > members;
    GetOutputSymbols(fst, &members);

    // for each unassigned label
    for (size_t label = 0; label < assigned_labels.size(); label++) {
      if (not assigned_labels[label]) {

        // reset the best hyp
        vector<VLabel> best_hyp(assigned_labels.size(), fst::SymbolTable::kNoSymbol);
        float max_member_score = Arc::Weight::Zero().Value();
        // for each possible member for this label
        for (typename std::set<typename Arc::Label>::const_iterator mem_it = members[label].begin();
            mem_it != members[label].end(); ++mem_it)
        {
          size_t member = *mem_it;

          // build a fst where (label, member) is selected
          // and some arcs are disabled according to mapper
//          Filter filter(label, member);
//          Mapper mapper(filter, dead_state);
//          fst::MapFst<Arc, Arc, Mapper> _fst(fst, mapper);

//          float score = -FstMass(_fst);

          // accumulate the score for this label
//          scores[label] = add_log(scores[label], score);
//          if (score > max_member_score) {
//            max_member_score = score;
//            copy(lhyp.begin(), lhyp.end(), best_hyp.begin());
//          }
        }


        // if the score for this label is better than for the rest, update
        if (scores[label] > max_score) {
          max_score = scores[label];
          copy(best_hyp.begin(), best_hyp.end(), hyp.begin());
        }
      } // if not assigned label
    }

  }
};

}

#endif /* openlat_INTERACTIVE_SEQUENCE_LABELING_H_ */
