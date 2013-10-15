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
 * lat-normalize.cpp
 *
 *  Created on: 15/03/2012
 *      Author: valabau
 */



#include <cmath>
#include <memory>
#include <sstream>

#include <fst/fstlib.h>
#include <fst/mutable-fst.h>
#include <openlat/compat.h>
#include <openlat/utils.h>
#include <openlat/normalize.h>
#include <openlat/iofilter.h>




using namespace std;
using namespace fst;
using namespace openlat;

typedef fst::VectorFst<fst::LogArc> LogVectorFst;

template <typename Arc>
void perform_normalization(const char* input, const char* output, float amscale, bool do_determinization) {
  {
    ifilter is(input);
    MutableFst <Arc> *fst = MutableFst <Arc> ::Read(is, FstReadOptions(input));
    Verify (*fst);
    if (amscale != 1.0) {
      ArcMap(fst, PowerMapper <Arc> (amscale));
    }
    if (do_determinization) {
      MutableFst <Arc> *ofst = new VectorFst<Arc>();
      DeterminizeAndNormalize(*fst, ofst);
      swap(fst, ofst);
      delete ofst;
    } else {
      Normalize(fst);
    }
    cerr << "fst is probabilistic = " << VerifyProbabilistic(*fst, 1e-4)
        << "\n";
    float entropy = Entropy(*fst);
    cerr << "ent = " << -entropy << "; ppl = " << exp(-entropy) << "\n";
    FstWriteOptions opts(output);
    ofilter os(output);
    fst->Write(os, opts);
    delete fst;
  }
}

int main(int argc, char *argv[]) {
  const string stdio_str("-");
  const char * input = stdio_str.c_str();
  const char *output = stdio_str.c_str();
  bool do_determinization = false, verbosity = false, show_help = false;
  bool is_std_fst = false;
  float amscale = 1.0;


  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "da:shv:")) != -1) {
    switch (c) {
      case 'd':
        do_determinization = true;
        break;
      case 's':
        is_std_fst = true;
        break;
      case 'a':
        amscale = convert_string<float>(optarg);
        break;
      case 'v':
        verbosity = convert_string<int>(optarg);
        break;
      case 'h':
        show_help = true;
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort();
        break;
    }
  }

  if (show_help) {
    cerr << "Usage: " << argv[0] << " [-d] [-a amscale] [-h] [-v int] [input-file] [output-file]\n";
    return 0;
  }

  if (optind < argc) input = argv[optind++];
  if (optind < argc) output = argv[optind++];

  if (verbosity && do_determinization)
    cerr << "Determinize\n";

  if (is_std_fst) {
    perform_normalization<StdArc>(input, output, amscale, do_determinization);
  }
  else {
    perform_normalization<LogArc>(input, output, amscale, do_determinization);
  }

  return EXIT_SUCCESS;
}

