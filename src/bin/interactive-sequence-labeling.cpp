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
#include <openlat/htk-compiler.h>
#include <openlat/iofilter.h>
#include <openlat/path-count.h>


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

  bool normalize = false;
  if (argc >= 7) normalize = (string(argv[6]) == "true");

  vector<const MutableFst<LogArc> *> fsts;
  vector<vector <VLabel> > refs;

  size_t n_lattices = 0;
  float total_paths = .0;
  float pruned_paths = .0;

  cerr << "Loading lattices with amscale " << amscale << " and pruning threshold " << pruning_threshold << " " << "with normalization = " << normalize << "\n";
  while (getline(refs_in, reference)) {
    file_in >> filename;


    //cerr << "Loading " << filename << " ...\n";
    VectorFst<LogArc> *fst = VectorFst<LogArc>::Read(filename);
    //ifilter is(filename.c_str());
    //MutableFst<LogArc> *fst = ReadHtkLogArc(is, filename);

    if (amscale != 1.0) ArcMap(fst, PowerMapper<LogArc>(amscale));

    n_lattices++;
    total_paths += PathCount(*fst); 
    if (pruning_threshold != 0) {
      PruneArcs(fst, pruning_threshold, 1e100);
      Connect(fst);
      if (normalize) Normalize(fst);
      float n_paths = PathCount(*fst); 
      if (n_paths < 1) {
        cerr << "Prunnig to heavy for lattice '" << filename << "'. Zero paths remain after the pruning.\n";
      }
      else {
        cerr << "lattice '" << filename << "' has " << n_paths << " paths\n";
      }
      pruned_paths += n_paths; 
    }

    vector<VLabel> ref;
    string_to_syms(reference, fst->OutputSymbols(), ref);
    //cerr << "with ref: " << syms_to_string(ref, fst->OutputSymbols()) << " ...\n";

    fsts.push_back(fst);
    refs.push_back(ref);
  }


  if (pruning_threshold != 0) {
    cerr << "Pruning to " << pruning_threshold << " results into " << (100.0*pruned_paths/total_paths) << "% of the paths\n";
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

  /* Pasive systems. */
  else if (method == "sequential-viterbi") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeViterbi<LogArc, LogConstraintFilter>, sort_pool_by_label>(fsts);
  }
  else if (method == "sequential-hamming") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter>, sort_pool_by_label>(fsts);
  }


  /* Active systems. Query strategies at structure level.  */
  else if (method == "active-structure-random") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label,
                             RandomSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-structure-entropy") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label,
                             EntropySampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-structure-least-confident") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-structure-error") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeSequentialExpectation<LogArc, LogConstraintFilter>, sort_pool_by_label,
                             ExpectedErrorSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }


  /* Active systems. Query strategies at label level.  */
  else if (method == "active-label-random") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             RandomLabelScorer<LogArc> >,       sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-entropy") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             EntropyLabelScorer<LogArc> >,       sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-least-confident") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             IdentityLabelScorer<LogArc> >,       sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-viterbi-error") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedShortestPathLabelScorer<LogArc> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-hamming-error") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedAccuracyLabelScorer<LogArc> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-mutual-information") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             MutualInformationLabelScorer<LogArc> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-num-changes") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestPathFunc<LogArc> > >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-hamming") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestHammingPathFunc<LogArc> > >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }

  /* Active systems. Query strategies at label level. More is less */
  else if (method == "active-label-entropy-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             EntropyLabelScorer<LogArc, false> >,       sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-least-confident-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             IdentityLabelScorer<LogArc, false> >,       sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-viterbi-error-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedShortestPathLabelScorer<LogArc, false> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-hamming-error-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedAccuracyLabelScorer<LogArc, false> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-mutual-information-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             MutualInformationLabelScorer<LogArc, false> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestPathFunc<LogArc>, false> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-hamming-m") {
    system = new LocalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestHammingPathFunc<LogArc>, false> >, sort_pool_by_score,
                             LeastConfidentSampleScorer<LogArc>, sort_sample_score_by_score>(fsts);
  }

  /* Active systems. Query strategies at label level with global pooling.  */
  else if (method == "active-label-random-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             RandomLabelScorer<LogArc> >,       sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-entropy-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             EntropyLabelScorer<LogArc> >,       sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-least-confident-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             IdentityLabelScorer<LogArc> >,       sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-viterbi-error-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedShortestPathLabelScorer<LogArc> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-hamming-error-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedAccuracyLabelScorer<LogArc> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-mutual-information-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             MutualInformationLabelScorer<LogArc> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestPathFunc<LogArc> > >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-hamming-global") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestHammingPathFunc<LogArc> > >, sort_pool_by_score>(fsts);
  }

  /* Active systems. Query strategies at label level with global pooling. More is less */
  else if (method == "active-label-entropy-global-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             EntropyLabelScorer<LogArc, false> >,       sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-least-confident-global-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             IdentityLabelScorer<LogArc, false> >,       sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-viterbi-error-global-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedShortestPathLabelScorer<LogArc, false> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-hamming-error-global-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedAccuracyLabelScorer<LogArc, false> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-global-mutual-information-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             MutualInformationLabelScorer<LogArc, false> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-global-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestPathFunc<LogArc>, false> >, sort_pool_by_score>(fsts);
  }
  else if (method == "active-label-num-changes-hamming-global-m") {
    system = new GlobalSystem<LogArc, LogConstraintFilter, RecomputeExpectedHamming<LogArc, LogConstraintFilter,
                             ConditionalExpectedNumChangesLabelScorer<LogArc, ShortestHammingPathFunc<LogArc>, false> >, sort_pool_by_score>(fsts);
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


