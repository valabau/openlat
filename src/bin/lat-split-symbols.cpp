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
#include <unistd.h>


#include <fst/fstlib.h>
#include <fst/mutable-fst.h>
#include <openlat/compat.h>
#include <openlat/utils.h>
#include <openlat/split-symbols.h>
#include <openlat/iofilter.h>
#include <fst/script/fst-class.h>



using namespace std;
using namespace fst;
using namespace fst::script;
using namespace openlat;

typedef fst::VectorFst<fst::LogArc> LogVectorFst;

template <typename Arc>
void process(const FstClass& _fst, const char *output, const string& separator, const string* space) {
  const Fst<Arc>& fst = *_fst.GetFst<Arc>();
  Verify(fst);

  fst::SymbolTable * const_syms = new fst::SymbolTable("const syms");
  const_syms->AddSymbol("<s>");
  const_syms->AddSymbol("</s>");
  const_syms->AddSymbol("<space>");
  const_syms->AddSymbol("<phrase>");
  const_syms->AddSymbol("<epsilon>");
  const_syms->AddSymbol("!NULL");

  VectorFst<Arc> ofst;
  SplitSymbols<Arc>(fst, &ofst, separator, space, const_syms);

  delete const_syms;

  FstWriteOptions opts(output);
  ofilter os(output);
  ofst.Write(os, opts);
}


int main(int argc, char *argv[]) {
  const string stdio_str("-");
  const char * input = stdio_str.c_str();
  const char *output = stdio_str.c_str();
  int verbosity = 0;
  bool show_help = false;

  string separator = "_";
  string *space = 0;

  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "s:a:v:h")) != -1) {
    switch (c) {
      case 's':
        separator = optarg;
        break;
      case 'a':
        space = new string(optarg);
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
    cerr << "Usage: " << argv[0] << " [-s separator_string ] [-a string_to_append_at_end ] [-h] [-v int] [input-file] [output-file]\n";
    return 0;
  }

  if (optind < argc) input = argv[optind++];
  if (optind < argc) output = argv[optind++];

  if (verbosity) {
    cerr << "separator = \"" << separator << "\"";
    if (space) cerr << " add \"" << *space << "\" at the end";
    cerr << "\n";
  }

  {
    ifilter is(input);
    FstClass *fst = FstClass::Read(is, string(input));

    const string& type = fst->ArcType();
    if (type == "log") {
      process<LogArc>(*fst, output, separator, space);
    }

    else if (type == "standard") {
      process<StdArc>(*fst, output, separator, space);
    }

    else {
      cerr << "Unsupported arc type '" << type << "'\n";
    }

    delete fst;
  }

  if (space) delete space;

  return EXIT_SUCCESS;
}

