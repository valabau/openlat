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



class System {
public:
  vector<LogVectorFst *> fsts;
  unordered_map<SampleLabel, PoolV> elems;
  vector< vector<bool> > assigned_labels;
  System(vector<LogVectorFst *> _fsts): fsts(_fsts) {
    assigned_labels.resize(fsts.size());
    for (size_t s = 0; s < fsts.size(); s++) {
      size_t n_labels = fsts[s]->InputSymbols()->NumSymbols();
      assigned_labels[s].resize(fsts[s]->InputSymbols()->NumSymbols());
      for (size_t l = 0; l < n_labels; l++) {
        elems.insert(make_pair<SampleLabel, PoolV>(SampleLabel(s, l), PoolV(0, 0)));
      }
    }

  }

  virtual Query askLabel() = 0;
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) = 0;
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) = 0;
  virtual void acceptLabel(const Query &query) {
    assigned_labels[query.label.sample][query.label.label] = true;

    LogVectorFst * fst = fsts[query.label.sample];
    LogQueryFilter filter(query.label.label, query.hyp);
    RmArc<LogArc, LogQueryFilter>(fst, RmArcOptions<LogArc, LogQueryFilter>(filter));
    assert_bt(fst->NumStates() > 0, "ERROR: empty lattice after accepting label\n");

  }
  virtual void update(const Query &query, vector<VLabel> &hyp) {
    if (query.label.label != kNoLabel) {
      acceptLabel(query);
    }

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

class LocalSystem: public System {
public:
  size_t next_sample;
  vector< vector<Pool::iterator>  > pool;
  LocalSystem(vector<LogVectorFst *>& _fsts): System(_fsts), next_sample(0), pool(fsts.size()) {
    for (Pool::iterator it = elems.begin(); it != elems.end(); ++it) {
      pool[it->first.sample].push_back(it);
    }
  }

  virtual Query askLabel() {
    Query query;
    query.label.sample = next_sample;
    next_sample = (next_sample + 1) % fsts.size();
    // pop min element from pool[query.sample]
    pop_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), less<Pool::iterator>());
    // return query
    query.label.label = pool[query.label.sample].back()->first.label;
    query.hyp         = pool[query.label.sample].back()->second.hyp;
    pool[query.label.sample].pop_back();
    return query;
  }
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    update(query, hyp);
    make_heap(pool[query.label.sample].begin(), pool[query.label.sample].end(), less<Pool::iterator>());
  }
  virtual ~LocalSystem() {};
};

class GlobalSystem: public System {
public:
  vector<Pool::iterator> pool;
  GlobalSystem(vector<LogVectorFst *>& _fsts): System(_fsts) {
    for (Pool::iterator it = elems.begin(); it != elems.end(); ++it) {
      pool.push_back(it);
    }
  }

  virtual Query askLabel() {
    Query query;
    // pop min element from pool
    pop_heap(pool.begin(), pool.end(), less<Pool::iterator>());
    // return query
    query.label.sample = pool.back()->first.sample;
    query.label.label = pool.back()->first.label;
    query.hyp         = pool.back()->second.hyp;
    pool.pop_back();
    return query;
  }
  virtual void fixLabel(const Query &query, vector<VLabel> &hyp) {
    update(query, hyp);
    make_heap(pool.begin(), pool.end(), less<Pool::iterator>());
  }
  virtual ~GlobalSystem() {};
};



void recompute_random_karyo(const LogVectorFst &fst, const vector<bool>&, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
  ShortestPath(fst, hyp, scores);
  for (size_t l = 0; l < scores.size(); l++) {
    scores[l] = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
  }

}

void recompute_greedy_karyo(const LogVectorFst &fst, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
  ShortestPath(fst, hyp, scores);
}

void recompute_l2r_karyo(const LogVectorFst &fst, const vector<bool>&, vector<VLabel> &hyp, vector<float> &scores) {
  ShortestPath(fst, hyp, scores);
  for (size_t l = 0; l < hyp.size(); l++) {
    scores[l] = -l;
  }
}

void recompute_global_karyo(const LogVectorFst &fst, LogVectorFst::StateId dead_state, const vector<bool>& assigned_labels,
                            vector<VLabel> &hyp, vector<float> &scores)
{
  scores.clear();
  scores.resize(assigned_labels.size(), -numeric_limits<float>::infinity());

  hyp.resize(assigned_labels.size());

//  cout << "assigned labels " << assigned_labels.to_string() << "\n";
//  cout << "assigned members " << assigned_members.to_string() << "\n";

  float max_score = -INFINITY;
  for (size_t label = 0; label < assigned_labels.size(); label++) {
    vector<VLabel> best_hyp(assigned_labels.size(), SymbolTable::kNoSymbol);

    if (not assigned_labels[label]) {
      float max_img_score = -INFINITY;
      for (size_t member = 0; member < fst.OutputSymbols()->NumSymbols(); member++) {
        LogQueryFilter filter(label, member);
        QueryMapper mapper(filter, dead_state);

        MapFst<LogArc, LogArc, QueryMapper> _fst(fst, mapper);
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


class RandomKaryoLocal: public LocalSystem {
public:
  RandomKaryoLocal(vector<LogVectorFst *>& _fsts): LocalSystem(_fsts) {};
  virtual ~RandomKaryoLocal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    recompute_greedy_karyo(*fsts[query.label.sample],
        assigned_labels[query.label.sample], hyp, scores);
  }
};

class RandomKaryoGlobal: public GlobalSystem {
public:
  RandomKaryoGlobal(vector<LogVectorFst *>& _fsts): GlobalSystem(_fsts) {};
  virtual ~RandomKaryoGlobal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    recompute_greedy_karyo(*fsts[query.label.sample],
        assigned_labels[query.label.sample], hyp, scores);
  }
};

class GreedyKaryoLocal: public LocalSystem {
public:
  GreedyKaryoLocal(vector<LogVectorFst *>& _fsts): LocalSystem(_fsts) {};
  virtual ~GreedyKaryoLocal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    recompute_greedy_karyo(*fsts[query.label.sample],
        assigned_labels[query.label.sample], hyp, scores);
  }
};

class GreedyKaryoGlobal: public GlobalSystem {
public:
  GreedyKaryoGlobal(vector<LogVectorFst *>& _fsts): GlobalSystem(_fsts) {};
  virtual ~GreedyKaryoGlobal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    recompute_greedy_karyo(*fsts[query.label.sample],
        assigned_labels[query.label.sample], hyp, scores);
  }
};

class L2RKaryoLocal: public LocalSystem {
public:
  L2RKaryoLocal(vector<LogVectorFst *>& _fsts): LocalSystem(_fsts) {};
  virtual ~L2RKaryoLocal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    recompute_l2r_karyo(*fsts[query.label.sample],
        assigned_labels[query.label.sample], hyp, scores);
  }
};

class GlobalKaryoLocal: public LocalSystem {
public:
  GlobalKaryoLocal(vector<LogVectorFst *>& _fsts): LocalSystem(_fsts) {};
  virtual ~GlobalKaryoLocal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    LogVectorFst *fst = fsts[query.label.sample];
    LogVectorFst::StateId dead_state = fst->AddState();
    recompute_global_karyo(*fst, dead_state,
        assigned_labels[query.label.sample], hyp, scores);
    fst->DeleteStates(vector<LogVectorFst::StateId>(1, dead_state));
  }
};

class GlobalKaryoGlobal: public GlobalSystem {
public:
  GlobalKaryoGlobal(vector<LogVectorFst *>& _fsts): GlobalSystem(_fsts) {};
  virtual ~GlobalKaryoGlobal() {};
  virtual void recomputeScores(const Query &query, vector<VLabel> &hyp, vector<float> &scores) {
    LogVectorFst *fst = fsts[query.label.sample];
    LogVectorFst::StateId dead_state = fst->AddState();
    recompute_global_karyo(*fst, dead_state,
        assigned_labels[query.label.sample], hyp, scores);
    fst->DeleteStates(vector<LogVectorFst::StateId>(1, dead_state));
  }
};


double diffclock(clock_t clock1,clock_t clock2) {
  double diffticks = clock1 - clock2;
  double diffms = (diffticks * 1000) / CLOCKS_PER_SEC;
  return diffms;
}

class Oracle {
public:
  const vector<vector<VLabel> >& refs;
  vector<size_t> n_errors;
  size_t c; // corrections
  size_t e_c; // label errors
  size_t e_k; // karyotype errors
  clock_t start_time;
  clock_t end_time;
  Oracle(const vector<vector<VLabel> >& _refs):
          refs(_refs), n_errors(refs.size()), c(0)
  {
    for (size_t i = 0; i < refs.size(); i++) n_errors[i] = refs[i].size();
    recomputeErrors();
  }

  void recomputeErrors() {
    e_c = 0;
    e_k = 0;
    for (size_t i = 0; i < n_errors.size(); i++) {
      if (n_errors[i] > 0) e_k++;
      e_c += n_errors[i];
    }
  }

  void printStats(size_t s) {
    recomputeErrors();
    end_time = clock();
    cout << s << " " << c << " " << e_c << " " << e_k << " " << diffclock(end_time, start_time)/1000 << endl;
  }

  void evaluate(System *system) {
    start_time = clock();
    for (size_t s = 0; s < refs.size(); s++) {
      Query query(s, SymbolTable::kNoSymbol, SymbolTable::kNoSymbol);
      vector<VLabel> hyp;
      system->fixLabel(query, hyp);
      n_errors[query.label.sample] = edit_distance(refs[query.label.sample], hyp);
    }
    printStats(0);
    size_t total_labels = 0;
    for (size_t i = 0; i < refs.size(); i++) total_labels += refs[i].size();
    for (size_t s = 0; s < total_labels; s++) {
      Query query = system->askLabel();
//      cout << "query: " << query.label.sample << " " << query.label.label << " " << query.hyp << " is " << refs[query.label.sample][query.label.label] << "\n";
      if (query.hyp != refs[query.label.sample][query.label.label]) {
        c++;
        query.hyp = refs[query.label.sample][query.label.label];
        vector<VLabel> hyp;
        system->fixLabel(query, hyp);
        n_errors[query.label.sample] = edit_distance(refs[query.label.sample], hyp);
      }
      else {
//        system->acceptLabel(query);
        vector<VLabel> hyp;
        system->fixLabel(query, hyp);
      }
      printStats(s + 1);
    }
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

  vector<LogVectorFst *> fsts;
  vector<vector <VLabel> > refs;

  while (getline(refs_in, reference)) {
    file_in >> filename;


    cout << "Loading " << filename << " ...\n";
    LogVectorFst *fst = LogVectorFst::Read(filename);

// XXX:   fst->scorer.amscale = amscale;

// XXX:   if (pruning_threshold != 0) {
//      LogVectorFst * fst_pruned = new LogVectorFst();
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
//  RandomKaryoLocal system(fsts);
  GlobalKaryoLocal system(fsts);
//  GlobalKaryoGlobal system(fsts);
//  L2RKaryoLocal system(fsts);

  Oracle oracle(refs);

  oracle.evaluate(&system);
  return EXIT_SUCCESS;

}

