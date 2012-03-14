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
 * htk-compiler.cpp
 *
 *  Created on: 23/02/2012
 *      Author: valabau
 */

#include "htk.h"
#include <openlat/htk-compiler.h>

using namespace fst;

namespace openlat {

template <typename Arc>
MutableFst<Arc>*        ReadHtk(std::istream &istrm, const std::string &source) {
  HtkContext context(istrm, source);

  int status = htk_parse(&context);

  if (status == 1) LOG(FATAL) << " Parsing '" << source << "' failed because of invalid input ";
  else if (status == 2) LOG(FATAL) << " Parsing '" << source << "' failed due to memory exhaustion ";

  return context.htk.CreateFst<Arc>();
}

fst::MutableFst<fst::LogArc>* ReadHtkLogArc(std::istream &istrm, const std::string &source) { return ReadHtk<fst::LogArc>(istrm, source); }
fst::MutableFst<fst::StdArc>* ReadHtkStdArc(std::istream &istrm, const std::string &source) { return ReadHtk<fst::StdArc>(istrm, source); }
fst::MutableFst<LogLinearArc>* ReadHtkLogLinearArc(std::istream &istrm, const std::string &source) { return ReadHtk<LogLinearArc>(istrm, source); }




}
