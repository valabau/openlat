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
#include <openlat/compat.h>
#include <openlat/utils.h>
#include <openlat/htk-compiler.h>
#include <openlat/iofilter.h>




using namespace std;
using namespace fst;
using namespace openlat;

extern int htk_debug;

int main(int argc, char *argv[]) {
  const string stdio_str("-");
  const char * input = stdio_str.c_str();
  const char *output = stdio_str.c_str();
  const char *epsilons_fn = 0;
  Wordlist epsilons;
  bool do_standard = false;

  if (argc >= 2 and argv[1][0] == '-' and argv[1][1] == 's') {
    argc--;
    argv++;
    do_standard = true;
  }

  if (argc >= 2) input = argv[1];
  if (argc >= 3) output = argv[2];
  if (argc >= 4) epsilons_fn = argv[3];

  
  if (epsilons_fn) {
    std::ifstream epsilons_s(epsilons_fn);
    std::string line;
    while (std::getline(epsilons_s, line)) {
      epsilons.insert(trim(line));
    }
  }
  

  if (do_standard) {
    ifilter is(input);
    MutableFst<StdArc> *fst = ReadHtkStdArc(is, input, epsilons);
    Verify(*fst);

    FstWriteOptions opts(output);
    ofilter os(output);
    fst->Write(os, opts);

    delete fst;
  }
  else {
    ifilter is(input);
    MutableFst<LogArc> *fst = ReadHtkLogArc(is, input, epsilons);
    Verify(*fst);

    FstWriteOptions opts(output);
    ofilter os(output);
    fst->Write(os, opts);

    delete fst;
  }

  return EXIT_SUCCESS;
}

