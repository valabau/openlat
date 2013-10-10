// compile.h

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Author: valabau@gmail.com (Vicent Alabau)
//
// \file
// Class to to compile a binary Fst from a lattice in THOT format.

#ifndef THOT_COMPILER_HPP_
#define THOT_COMPILER_HPP_

#include <string>
#include <vector>

#include <iostream>
#include <fstream>

#include <openlat/htk-compiler.h>

namespace openlat {

fst::MutableFst<fst::LogArc>*   ReadThotLogArc(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols = Wordlist(), int debug_level = -1);
fst::MutableFst<fst::StdArc>*   ReadThotStdArc(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols = Wordlist(), int debug_level = -1);
fst::MutableFst<LogLinearArc>*  ReadThotLogLinearArc(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols = Wordlist(), int debug_level = -1);

Lattice<fst::LogArc>*  ReadThotLogLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols = Wordlist(), int debug_level = -1);
Lattice<fst::StdArc>*  ReadThotStdLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols = Wordlist(), int debug_level = -1);
Lattice<LogLinearArc>* ReadThotLogLinearLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols = Wordlist(), int debug_level = -1);

}  // namespace fst

#endif  //THOT_COMPILER_HPP_
