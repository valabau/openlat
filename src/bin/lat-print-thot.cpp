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

#include <fst/fstlib.h>
#include <fst/script/fst-class.h>
#include <openlat/compat.h>
#include <openlat/utils.h>
#include <openlat/thot-printer.h>
#include <openlat/iofilter.h>



using namespace std;
using namespace fst;
using namespace openlat;

template<class Arc>
void execute(const fst::Fst<Arc> &fst, const char *output, const char *start_symbol, const char *end_symbol) {
  Verify(fst);

  ofilter os(output);

  if (start_symbol or end_symbol) {
    MutableFst<Arc> *mfst = new VectorFst<Arc>(fst);

//    AddStartAndEndSymbols(mfst, start_symbol, end_symbol);
    PrintThot(*mfst, os);

    delete mfst;
  }
  else {
    PrintThot(fst, os);
  }
}

int main(int argc, char *argv[]) {
  const string stdio_str("-");
  const char * input = stdio_str.c_str();
  const char *output = stdio_str.c_str();
  char *start_symbol = 0;
  char *end_symbol = 0;

  if (argc >= 4 and argv[1][0] == '-' and argv[1][1] == 's') {
    start_symbol = argv[2];
    end_symbol   = argv[3];
    argc -= 3;
    argv += 3;
  }

  if (argc >= 2)  input = argv[1];
  if (argc >= 3) output = argv[2];

  script::FstClass *fst = script::FstClass::Read((strncmp(input, "-", 1) == 0)?"":input);

  if (not fst) return EXIT_FAILURE;

  if (fst->ArcType() == "standard") {
    execute(*fst->GetFst<StdArc>(), output, start_symbol, end_symbol);
  }
  else if (fst->ArcType() == "log") {
    execute(*fst->GetFst<LogArc>(), output, start_symbol, end_symbol);
  }

  delete fst;

  return EXIT_SUCCESS;
}

