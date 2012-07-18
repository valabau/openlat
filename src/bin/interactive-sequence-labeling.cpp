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
#include <openlat/normalize.h>


#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <ctime>

using namespace fst;
using namespace openlat;


int main (int argc, char *argv[]) {
  //INIT_TRACE();

  cerr.precision(3);

  assert_bt(argc >= 4, "Invalid number of parameters");
  ifstream file_in(argv[2]);
  ifstream refs_in(argv[3]);
  string reference;
  string filename;

  float amscale = 1.0;
  if (argc >= 5) amscale = convert_string<float>(argv[4]);

  float pruning_threshold = .0;
  if (argc >= 6) pruning_threshold = convert_string<float>(argv[5]);

  vector<const MutableFst<LogArc> *> fsts;
  vector<vector <VLabel> > refs;

  cerr << "Loading lattices\n";
  while (getline(refs_in, reference)) {
    file_in >> filename;


    //cerr << "Loading " << filename << " ...\n";
    VectorFst<LogArc> *fst = VectorFst<LogArc>::Read(filename);

    if (amscale != 1.0) ArcMap(fst, PowerMapper<LogArc>(amscale));

// XXX:   if (pruning_threshold != 0) {
//      VectorFst<Arc> * fst_pruned = new VectorFst<Arc>();
//      fst->copyPrunePosterior(*fst_pruned, pruning_threshold);
//      swap(fst, fst_pruned);
//      delete fst_pruned;
//    }

    vector<VLabel> ref;
    string_to_syms(reference, fst->OutputSymbols(), ref);
    //cerr << "with ref: " << syms_to_string(ref, fst->OutputSymbols()) << " ...\n";

    fsts.push_back(fst);
    refs.push_back(ref);
  }

  assert_bt(fsts.size() > 0, "No fsts exist!!!");

  cerr << "Finish loading lattices\n";

  srand(time(NULL));

  string method(argv[1]);

  ISystem *system = 0;

  if (method == "random") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeRandom<LogArc, LogConstraintFilter> >(fsts);
  }
  else if (method == "sequential-path") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequential<LogArc, LogConstraintFilter>, sort_pool_by_label>(fsts);
  }
  else if (method == "sequential-symbol") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label>(fsts);
  }
  else if (method == "sequential-viterbi") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeViterbi<LogArc, LogConstraintFilter>, sort_pool_by_label>(fsts);
  }
  else if (method == "sequential-hamming") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter>, sort_pool_by_label>(fsts);
  }
  else if (method == "sorted-viterbi") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeViterbi<LogArc, LogConstraintFilter>, sort_pool_by_label, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-viterbi-hamming") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeViterbi<LogArc, LogConstraintFilter, ExpectedHammingLabelScorer<LogArc> >, sort_pool_by_label, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-sequential-symbol") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-sequential-hamming") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter>,       sort_pool_by_label, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-hamming") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter>,       sort_pool_by_score, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-hamming-num-changes") {
    system = new LocalSystem<LogArc, LogConstraintFilter,
                             RecomputeExpectedHamming<LogArc, LogConstraintFilter, ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestPathFunc<LogArc> > >,
                             sort_pool_by_score, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-hamming-mutual-information") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter, MutualInformationLabelScorer<LogArc> >, sort_pool_by_score, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-hamming-expected-accuracy") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter, ConditionalExpectedAccuracyLabelScorer<LogArc> >, sort_pool_by_score, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-viterbi-expected-accuracy") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter, ConditionalExpectedAccuracyLabelScorer<LogArc> >, sort_pool_by_score, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-sequential-symbol-posterior") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_score, EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-sequential-mutual-information") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label, MaxMutualInformationSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-sequential-expected-error") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label, ExpectedErrorSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "sorted-updated-sequential-symbol") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label, EntropySampleScorer<LogArc>, sort_sample_score_by_score, true>(fsts);
  }
  else if (method == "sorted-updated-sequential-mutual-information") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label, MaxMutualInformationSampleScorer<LogArc>, sort_sample_score_by_score, true>(fsts);
  }
  else if (method == "sorted-updated-sequential-expected-error") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label, ExpectedErrorSampleScorer<LogArc>, sort_sample_score_by_score, true>(fsts);
  }
  else if (method == "reverse-sequential-path") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequential<LogArc, LogConstraintFilter>, sort_pool_by_label_reverse>(fsts);
  }
  else if (method == "reverse-sequential-symbol") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter, true>, sort_pool_by_label_reverse>(fsts);
  }
  else if (method == "sequential-greedy") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter, false, false>, sort_pool_by_label>(fsts);
  }
  else if (method == "reverse-sequential-greedy") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter, true, false>, sort_pool_by_label_reverse>(fsts);
  }
  else if (method == "active-path") {
    system = new TraversalLocalSystem<LogArc, LogConstraintFilter, RecomputeSequential<LogArc, LogConstraintFilter> >(fsts);
  }

  assert_bt(system != 0, "Invalid method name");

  const SymbolTable *isymbols = (*fsts.begin())->InputSymbols();
  const SymbolTable *osymbols = (*fsts.begin())->OutputSymbols();
  Oracle oracle(refs, isymbols, osymbols);

  oracle.evaluate(system);

  cout << oracle.c << "\n";

  delete system;
  for (size_t i = 0; i < fsts.size(); i++) delete fsts[i];
  return EXIT_SUCCESS;

}


