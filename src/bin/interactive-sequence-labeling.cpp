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


#include <fst/fstlib.h>
#include <openlat/interactive-sequence-labeling.h>


#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <ctime>

using namespace fst;
using namespace openlat;


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


