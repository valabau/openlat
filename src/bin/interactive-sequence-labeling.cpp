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


int main (int argc, char *argv[]) {
  //INIT_TRACE();

  cout.precision(3);

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

  cerr << "Loading lattices\n";
  while (getline(refs_in, reference)) {
    file_in >> filename;


    //cerr << "Loading " << filename << " ...\n";
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
    //cerr << "with ref: " << syms_to_string(ref, fst->OutputSymbols()) << " ...\n";

    fsts.push_back(fst);
    refs.push_back(ref);
  }

  cerr << "Finish loading lattices\n";

  srand(time(NULL));
//  LocalSystem<LogArc, LogQueryFilter, RecomputeGreedy<LogArc, LogQueryFilter> > system(fsts);
  LocalSystem<LogArc, LogQueryFilter, RecomputeSequential<LogArc, LogQueryFilter>, sort_pool_by_label> system(fsts);
//  LocalSystem<LogArc, LogQueryFilter, RecomputeExpectation<LogArc, LogQueryFilter> > system(fsts);
//  GlobalKaryoGlobal system(fsts);
//  LocalSystem<LogArc, LogQueryFilter, RecomputeSequential<LogArc, LogQueryFilter> > system(fsts);

  Oracle oracle(refs);

  oracle.evaluate(&system);

  for (size_t i = 0; i < fsts.size(); i++) delete fsts[i];
  return EXIT_SUCCESS;

}


