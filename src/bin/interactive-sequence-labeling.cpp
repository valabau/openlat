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
 * interactive-sequence-labeling.cpp
 *
 *  Created on: 17/02/2012
 *      Author: valabau
 */




#include <fst/interval-set.h>
#include <fst/fstlib.h>
#include <openlat/compat.h>
#include <openlat/rmarc.h>
#include <openlat/utils.h>
#include <openlat/query.h>


#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <ctime>

using namespace fst;
using namespace openlat;


typedef struct _sPoolV {
  VLabel hyp;
  float score;
  _sPoolV(VLabel _hyp, float _score): hyp(_hyp), score(_score) {}
  bool operator<(const _sPoolV &other) const { return score < other.score; }
} PoolV;

typedef unordered_map<SampleLabel, PoolV> Pool;

namespace std {
  template<> struct less<Pool::iterator> {
    size_t operator()(const Pool::iterator& __x, const Pool::iterator& __y) const {
      return __x->second < __y->second;
    }
  };
}

/**
 * @param clock1 first clock
 * @param clock2 second clock
 * @return return the lapse of time in ms between clock1 and clock2
 */
double diffclock(clock_t clock1, clock_t clock2) {
  double diffticks = clock1 - clock2;
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
      Query query(s, SymbolTable::kNoSymbol, SymbolTable::kNoSymbol);
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
  vector<VectorFst<Arc> *> fsts; //< store all fsts samples
  unordered_map<SampleLabel, PoolV> elems; //< pool of hypotheses
  vector< vector<bool> > assigned_labels; //< store already assigned labels for each sample

  System(vector<VectorFst<Arc> *> _fsts): fsts(_fsts) {
    assigned_labels.resize(fsts.size());
    for (size_t s = 0; s < fsts.size(); s++) {
      // XXX: it should be checked if all fst comply with the sequence labeling restrictions
      // obtain the number of labels for this sample
      // XXX: do not rely on the number of input symbols. make an explicit check instead
      size_t n_labels = fsts[s]->InputSymbols()->NumSymbols();
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
    VectorFst<Arc> * fst = fsts[query.label.sample];
    Filter filter(query.label.label, query.hyp);
    RmArc<Arc, Filter>(fst, RmArcOptions<Arc, Filter>(filter));
    assert_bt(fst->NumStates() > 0, "ERROR: empty lattice after accepting label\n");
  }

  /** update the system with the query and return the updated hypothesis */
  virtual void update(const Query &query, vector<VLabel> &hyp) {
    // accept if the query is known
    if (query.label.label != kNoLabel) {
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
  LocalSystem(vector<VectorFst<Arc> *>& _fsts): Base(_fsts), next_sample(0), pool(Base::fsts.size()) {
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
    pop_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), less<Pool::iterator>());

    // return query and remove it from the queue
    query.label.label = pool[query.label.sample].back()->first.label;
    query.hyp         = pool[query.label.sample].back()->second.hyp;
    pool[query.label.sample].pop_back();
    return query;
  }

  /** Assign query and update re-sort the queue jut for the current sample */
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    Base::update(query, hyp);
    make_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), less<Pool::iterator>());
  }

  /** Recompute scores using the Recompute template argument */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    VectorFst<Arc> *fst = Base::fsts[query.label.sample];
    // create a dummy dead state for the mappers to disable an arc
    typename VectorFst<Arc>::StateId dead_state = fst->AddState();
    Recompute()(*fst, dead_state,
        Base::assigned_labels[query.label.sample], hyp, scores);
    // delete the dummy state
    fst->DeleteStates(vector<typename VectorFst<Arc>::StateId>(1, dead_state));
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
  GlobalSystem(vector<VectorFst<Arc> *>& _fsts): Base(_fsts) {
    for (Pool::iterator it = Base::elems.begin(); it != Base::elems.end(); ++it) {
      pool.push_back(it);
    }
  }

  /** Ask next label */
  virtual Query askLabel() {
    Query query;
    // pop min element from the global pool
    pop_heap(pool.begin(), pool.end(), less<Pool::iterator>());
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
    make_heap(pool.begin(), pool.end(), less<Pool::iterator>());
  }

  /** Recompute scores using the Recompute template argument */
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    VectorFst<Arc> *fst = Base::fsts[query.label.sample];
    // create a dummy dead state for the mappers to disable an arc
    typename VectorFst<Arc>::StateId dead_state = fst->AddState();
    Recompute()(*fst, dead_state,
        Base::assigned_labels[query.label.sample], hyp, scores);
    // delete the dummy state
    fst->DeleteStates(vector<typename VectorFst<Arc>::StateId>(1, dead_state));
  }
  virtual ~GlobalSystem() {};

protected:
  vector<Pool::iterator> pool; //< a global pool of queries
};

template <class Arc, class Filter>
struct RecomputeRandom {
  void operator()(const VectorFst<Arc> &fst, typename VectorFst<Arc>::StateId, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    ShortestPath(fst, hyp, scores);
    for (size_t l = 0; l < scores.size(); l++) {
      scores[l] = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
    }
  }
};

template <class Arc, class Filter>
struct RecomputeGreedy {
  void operator()(const VectorFst<Arc> &fst, typename VectorFst<Arc>::StateId, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    ShortestPath(fst, hyp, scores);
  }
};

template <class Arc, class Filter>
struct RecomputeSequential {
  void operator()(const VectorFst<Arc> &fst, typename VectorFst<Arc>::StateId, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
    ShortestPath(fst, hyp, scores);
    for (size_t l = 0; l < hyp.size(); l++) {
      scores[l] = -l;
    }
  }
};

template <class Arc, class Filter>
struct RecomputeExpectation {
  void operator()(const VectorFst<Arc> &fst, typename VectorFst<Arc>::StateId dead_state, const vector<bool>& assigned_labels, vector<VLabel> &hyp, vector<float> &scores) {
    typedef FilterMapper<Arc, Filter> Mapper;

    scores.clear();
    scores.resize(assigned_labels.size(), numeric_limits<float>::infinity());

    hyp.resize(assigned_labels.size());

  //  cout << "assigned labels " << assigned_labels.to_string() << "\n";

    float max_score = -INFINITY;
    for (size_t label = 0; label < assigned_labels.size(); label++) {
      vector<VLabel> best_hyp(assigned_labels.size(), SymbolTable::kNoSymbol);

      if (not assigned_labels[label]) {
        float max_img_score = -INFINITY;
        for (size_t member = 0; member < fst.OutputSymbols()->NumSymbols(); member++) {
          Filter filter(label, member);
          Mapper mapper(filter, dead_state);

          MapFst<Arc, Arc, Mapper> _fst(fst, mapper);
          //XXX: is this really necessary? -> mapper.SetDeadState(_fst.AddState());

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
          if (score > max_img_score) {
            max_img_score = score;
            copy(lhyp.begin(), lhyp.end(), best_hyp.begin());
          }
        }

  //      cout << "max score l" << label << ": " << scores[label] << " ";
  //      cout << sentence_to_string(best_hyp);
  //      cout << "\n";
        if (scores[label] > max_score) {
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



void string_to_syms(const string &str, const SymbolTable *symtab, vector<VLabel> &syms) {
  if (str == "") return;

  vector<string> words_str;
  tokenize<string>(str, words_str, " \t\n");

  syms.clear();
  syms.reserve(words_str.size());
  for (size_t i = 0; i < words_str.size(); i++) {
    VLabel sym = symtab->Find(words_str[i]);
    assert_bt(sym != SymbolTable::kNoSymbol, "Invalid symbol '" << words_str[i] << "'\n");
    syms.push_back(sym);
  }
}


string syms_to_string(vector<VLabel> &syms, const SymbolTable *symtab) {
  stringstream ss;

  if (syms.empty()) return ss.str();

  ss << symtab->Find(syms[0]);
  for (size_t i = 1; i < syms.size(); i++) {
    ss << " " << symtab->Find(syms[i]);
  }

  return ss.str();
}


int main (int argc, char *argv[]) {
  //INIT_TRACE();

  assert_bt(argc >= 3, "Invalid number of parameters");
  ifstream file_in(argv[1]);
  ifstream refs_in(argv[2]);
  string reference;
  string filename;

  float amscale = 1.0;
  if (argc >= 4) amscale = convert_string<float>(argv[3]);

  float pruning_threshold = .0;
  if (argc >= 5) pruning_threshold = convert_string<float>(argv[4]);

  vector<VectorFst<LogArc> *> fsts;
  vector<vector <VLabel> > refs;

  while (getline(refs_in, reference)) {
    file_in >> filename;


    cout << "Loading " << filename << " ...\n";
    VectorFst<LogArc> *fst = VectorFst<LogArc>::Read(filename);

// XXX:   fst->scorer.amscale = amscale;

// XXX:   if (pruning_threshold != 0) {
//      VectorFst<Arc> * fst_pruned = new VectorFst<Arc>();
//      fst->copyPrunePosterior(*fst_pruned, pruning_threshold);
//      swap(fst, fst_pruned);
//      delete fst_pruned;
//    }

    vector<VLabel> ref;
    string_to_syms(reference, fst->OutputSymbols(), ref);
    cout << "with ref: " << syms_to_string(ref, fst->OutputSymbols()) << " ...\n";

    fsts.push_back(fst);
    refs.push_back(ref);
  }

  srand(time(NULL));
//  LocalSystem<LogArc, LogQueryFilter, RecomputeGreedy<LogArc, LogQueryFilter> > system(fsts);
  LocalSystem<LogArc, LogQueryFilter, RecomputeExpectation<LogArc, LogQueryFilter> > system(fsts);
//  GlobalKaryoGlobal system(fsts);
//  LocalSystem<LogArc, LogQueryFilter, RecomputeSequential<LogArc, LogQueryFilter> > system(fsts);

  Oracle oracle(refs);

  oracle.evaluate(&system);
  return EXIT_SUCCESS;

}


