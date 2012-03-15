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

int main(int argc, char *argv[]) {
  const string stdio_str("-");
  const char * input = stdio_str.c_str();
  const char *output = stdio_str.c_str();

  if (argc >= 2) input = argv[1];
  if (argc >= 3) output = argv[2];
  {
    ifilter is(input);
    MutableFst<LogArc> *fst = MutableFst<LogArc>::Read(is, FstReadOptions(input));
    Verify(*fst);

    Normalize(fst);

    FstWriteOptions opts(output);
    ofilter os(output);
    fst->Write(os, opts);

    delete fst;
  }

  return EXIT_SUCCESS;
}

