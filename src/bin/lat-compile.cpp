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
 * lat-compile.cpp
 *
 *  Created on: 17/02/2012
 *      Author: valabau
 */


#include <cmath>
#include <memory>
#include <sstream>
#include <unistd.h>

#include <fst/fstlib.h>
#include <openlat/compat.h>
#include <openlat/utils.h>
#include <openlat/htk-compiler.h>
#include <openlat/thot-compiler.h>
#include <openlat/iofilter.h>




using namespace std;
using namespace fst;
using namespace openlat;

int main(int argc, char *argv[]) {
  const string stdio_str("-");
  const char * input = stdio_str.c_str();
  const char *output = stdio_str.c_str();
  const char *epsilons_fn = 0;
  Wordlist epsilons;
  bool do_standard = false;
  bool is_thot_wgp = false;
  bool show_help = false;
  int verbosity = 0;


  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "sthv:e:")) != -1) {
    switch (c) {
      case 's':
        do_standard = true;
        break;
      case 't':
        is_thot_wgp = true;
        break;
      case 'e':
        epsilons_fn = optarg;
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
    cerr << "Usage: " << argv[0] << " [-s] [-t] [-h] [-v int] [-e filename] [input-file] [output-file]\n";
    return 0;
  }

  if (optind < argc) input = argv[optind++];
  if (optind < argc) output = argv[optind++];

  
  if (epsilons_fn) {
    std::ifstream epsilons_s(epsilons_fn);
    std::string line;
    if (verbosity > 0) cerr << "epsilons:";
    while (std::getline(epsilons_s, line)) {
      if (verbosity > 0) cerr << " " << trim(line);
      epsilons.insert(trim(line));
    }
    if (verbosity > 0) cerr << "\n";
  }
  
  if (is_thot_wgp) cerr << "Reading thot wgp from '" << input << "'\n";
  else cerr << "Reading htk lattice from '" << input << "'\n";

  if (do_standard) {
    ifilter is(input);
    MutableFst<StdArc> *fst = 0;
    if (is_thot_wgp) {
#if HAVE_THOT
      fst = ReadThotStdArc(is, input, epsilons, verbosity);
#else
      cerr << "ERROR: openlat has not been compiled with Thot support.\n";
#endif
    }
    else {
      fst = ReadHtkStdArc(is, input, epsilons, verbosity);
    }
    Verify(*fst);

    FstWriteOptions opts(output);
    ofilter os(output);
    fst->Write(os, opts);

    delete fst;
  }
  else {
    ifilter is(input);
    MutableFst<LogArc> *fst = 0;
    if (is_thot_wgp) {
#if HAVE_THOT
      fst = ReadThotLogArc(is, input, epsilons, verbosity);
#else
      cerr << "ERROR: openlat has not been compiled with Thot support.\n";
#endif
    }
    else {
      fst = ReadHtkLogArc(is, input, epsilons, verbosity);
    }
    Verify(*fst);

    FstWriteOptions opts(output);
    ofilter os(output);
    fst->Write(os, opts);

    delete fst;
  }

  return EXIT_SUCCESS;
}

