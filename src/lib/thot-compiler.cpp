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
 * thot-compiler.cpp
 *
 *  Created on: 23/02/2012
 *      Author: valabau
 */

#include "thot/thot.h"
#include <openlat/thot-compiler.h>
#include <set>
#include <fstream>

namespace openlat {

namespace thot {
class ThotParser {
  ThotContext * context_;
  int debug_level_;
  size_t n_arcs_;
public:
  ThotParser(ThotContext *context): context_(context), debug_level_(0), n_arcs_(0) { }
  void set_debug_level(int debug_level) { debug_level_ = debug_level;  }
  int parse() {
    string line, superfinal = "***SUPERFINALSTATE***";
    getline(context_->is, line); // ignore header

    set<string> final_nodes;
    { // read nodes
      getline(context_->is, line);
      std::istringstream iss(line);
      string node_id;
      while (iss >> node_id) final_nodes.insert(node_id);
    }

    // read arcs

    while (getline(context_->is, line)) {
      std::istringstream iss(line);
      string start_node_id, end_node_id, words;
      float score;
      iss >> start_node_id >> end_node_id >> score;
      getline(iss, words);

      ThotLink *link = context_->thot.getNewArc();
      link->start = context_->thot.getStateId(start_node_id);
      if (final_nodes.find(end_node_id) != final_nodes.end()) end_node_id = superfinal;
      link->end   = context_->thot.getStateId(end_node_id);
      link->score = score;
      context_->thot.assignWords(link, trim(words));

      if (debug_level_ > 0) {
        cerr << n_arcs_ << " " << link->start << " " << link->end << " " << link->score << " '" << context_->thot.wsyms->Find(link->words) << "'\n";
      }
      n_arcs_++;
    }

    return 0;
  }
};
}

template <typename Arc>
MutableFst<Arc>* ReadThot(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level = -1) {
  thot::ThotContext context(istrm, source, epsilon_symbols);
  thot::ThotParser parser(&context);
  if (debug_level != -1) parser.set_debug_level(debug_level);

  int status = parser.parse();

  if (status == 1) LOG(FATAL) << " Parsing '" << source << "' failed because of invalid input ";
  else if (status == 2) LOG(FATAL) << " Parsing '" << source << "' failed due to memory exhaustion ";

  return context.thot.CreateFst<Arc>();
}

fst::MutableFst<fst::LogArc>* ReadThotLogArc(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level) { return ReadThot<fst::LogArc>(istrm, source, epsilon_symbols, debug_level); }
fst::MutableFst<fst::StdArc>* ReadThotStdArc(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level) { return ReadThot<fst::StdArc>(istrm, source, epsilon_symbols, debug_level); }
fst::MutableFst<LogLinearArc>* ReadThotLogLinearArc(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level) { return ReadThot<LogLinearArc>(istrm, source, epsilon_symbols, debug_level); }

template <typename Arc>
Lattice<Arc>* ReadThotLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level = -1) {
  thot::ThotContext context(istrm, source, epsilon_symbols);
  thot::ThotParser parser(&context);
  if (debug_level != -1) parser.set_debug_level(debug_level);

  int status = parser.parse();

  if (status == 1) LOG(FATAL) << " Parsing '" << source << "' failed because of invalid input ";
  else if (status == 2) LOG(FATAL) << " Parsing '" << source << "' failed due to memory exhaustion ";

  return context.thot.CreateLattice<Arc>();
}


Lattice<fst::LogArc>* ReadThotLogLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level) { return ReadThotLattice<fst::LogArc>(istrm, source, epsilon_symbols, debug_level); }
Lattice<fst::StdArc>* ReadThotStdLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level) { return ReadThotLattice<fst::StdArc>(istrm, source, epsilon_symbols, debug_level); }
Lattice<LogLinearArc>* ReadThotLogLinearLattice(std::istream &istrm, const std::string &source, const Wordlist& epsilon_symbols, int debug_level) { return ReadThotLattice<LogLinearArc>(istrm, source, epsilon_symbols, debug_level); }

}
