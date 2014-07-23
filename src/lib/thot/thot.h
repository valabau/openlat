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
 * thot.h
 *
 *  Created on: 23/02/2012
 *      Author: valabau
 */

#ifndef openlat_THOT_H_
#define openlat_THOT_H_

#include <string>
#include <vector>
#include <deque>
#include <unordered_map>

#include <fst/fst.h>
#include <fst/arc.h>
#include <fst/symbol-table.h>
#include <openlat/utils.h>
#include <openlat/vector-weight.h>
#include <openlat/thot-compiler.h>
#include <cmath>

using namespace std;
using namespace fst;


namespace openlat { namespace thot {

class ThotLink {
public:
  string name;     // Link name
  int start; // Start node number (of the link)
  int end;   // End node number (of the link)
  int words;   // Output symbol
  float score;
  vector<float> features; // Features

  ThotLink(): start(kNoStateId), end(kNoStateId), words(kNoLabel) {};
};

class ThotNode {
public:
  string name;     // Node name

  ThotLink link;    // auxiliar link
  ThotNode() {};
};

class ThotWeights {
public:
  typedef enum {XSCORE} feature_t;

  struct feature_hash {
    size_t operator()(const feature_t& __x) const {
      return std::hash<int>()(static_cast<int>(__x));
    }
  };

  ThotWeights(): amscale(1.0) {}

  vector<float> weights; // scales for features
  float  amscale;  // Acoustic scale for posterior computation
  unordered_map<feature_t, size_t, feature_hash> feature_pos; // position of the feature

  size_t getFeature(feature_t feat_name) {
    unordered_map<feature_t, size_t, feature_hash>::const_iterator it = feature_pos.find(feat_name);
    if (it == feature_pos.end()) {
      it = feature_pos.insert(make_pair(feat_name, feature_pos.size())).first;
    }
    return it->second;
  }

  std::vector<std::string> getFeatureNames() {
    unordered_map<feature_t, size_t, feature_hash>::const_iterator it;
    std::string strnames[] = {"score", "xscore" };

    std::vector<std::string> names; 
    for (it = feature_pos.begin(); it != feature_pos.end(); ++it) {
      if (it->second >= names.size()) names.resize(it->second + 1);
      if (it->first < XSCORE) names[it->second] = strnames[it->first];
      else names[it->second] = strnames[XSCORE] + to_string(it->first - XSCORE + 1);
    }
    return names;
  }

  std::vector<float> getWeights() {
    return weights;
  }

  size_t findFeature(feature_t feat_name) const {
    unordered_map<feature_t, size_t, feature_hash>::const_iterator it = feature_pos.find(feat_name);
    return (it == feature_pos.end()) ? static_cast<size_t>(-1): it->second;
  }

};

template <typename Weight>
struct compute_weight {
  inline Weight operator()(const ThotLink &link, const ThotWeights &weights) {
    float w = link.score;

    const size_t m = std::min(link.features.size(), weights.weights.size());
    if (m > 0) {
      size_t i;
      w = 0;
      for (i = 0; i < m; i++) {
        w += link.features[i] * weights.weights[i];
      }
      for (; i < link.features.size(); i++) {
        w += link.features[i];
      }
    }
    w /= -weights.amscale;
    return Weight(w);
  }
};

template <>
struct compute_weight<LogLinearWeight> {
  inline LogLinearWeight operator()(const ThotLink &link, const ThotWeights &weights) {
    if (link.features.empty()) {
      return LogLinearWeight::One();
    }
    std::vector<float> _w(link.features);
    for (size_t i = 0; i < _w.size(); i++) {
      _w[i] = -_w[i]; 
    }

    LogLinearWeight::Weight computed = compute_weight<LogLinearWeight::Weight>()(link, weights);
    LogLinearWeight w = LogLinearWeight(computed.Value(), _w);
    w.Resize(weights.feature_pos.size());
    return w;
  }
};



class ThotLattice {
public:

  string version;   // Lattice specification adhered to
  string utterance; // Utterance identifier
  float  base;      // LogBase for Likelihoods (0.0 not logs,default base e)
  string source;    // filename


  SymbolTable *ssyms;
  SymbolTable *wsyms;

  ThotWeights weights;
  vector<ThotLink> links;
  vector<ThotNode> nodes;

  string start_name;
  string end_name;

  const Wordlist& epsilon_symbols;

  ThotLattice(const Wordlist& _epsilon_symbols = Wordlist()): base(M_E),
                ssyms(new SymbolTable("states")), wsyms(new SymbolTable("words")),
                epsilon_symbols(_epsilon_symbols)
  {
    // add epsilon symbols (!NULL in thot)
    wsyms->AddSymbol("!NULL", 0);
    for (Wordlist::const_iterator it = _epsilon_symbols.begin(); it != _epsilon_symbols.end(); ++it) {
      wsyms->AddSymbol(*it, 0);
    }
  };

  ~ThotLattice() {
    delete ssyms; delete wsyms;
  };

  float getScore(const vector<float> scores) {
    float sum = .0;
    if (scores.size() != weights.weights.size()) LOG(FATAL) << "Score size and weigth size do not match\n";
    for (size_t i = 0; i < scores.size(); i++) {
      sum += scores[i] * weights.weights[i];
    }
    return sum / weights.amscale;
  }

  ThotLink *getNewArc() {
    size_t id = links.size();
    links.resize(id + 1);
    return &links[id];
  }

  ThotNode *getState(string state_name) {
    return &nodes[getStateId(state_name)];
  }

  LogLinearFst::StateId getStateId(string state_name) {
    size_t id = ssyms->AddSymbol(state_name);
    if (nodes.size() < id + 1) nodes.resize(id + 1);
    nodes[id].name = state_name;
    nodes[id].link.name = state_name + "_" + "aux";
    return id;
  }

  void assignWords(ThotLink *l, const string& output) {
    size_t id = wsyms->AddSymbol(output);
    if (output.empty() or epsilon_symbols.find(output) != epsilon_symbols.end()) {
      id = 0;
    }
    l->words = id;
  }

  int checkStartState() {
    int start_node = -1;
    if (start_name.empty()) {
      vector<size_t> n_in_arcs(nodes.size(), 0);
      for (vector<ThotLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
        n_in_arcs[link->end]++;
      }
      for (size_t n = 0; n < nodes.size(); n++) {
        if (n_in_arcs[n] == 0) {
          if (start_node == -1) start_node = n;
          else {
            LOG(FAIL) << " invalid lattice '" << source << "'; it has more than one initial state ";
          }
        }
      }
      if (start_node != -1) start_name = nodes[start_node].name;
    }
    else {
      start_node = ssyms->Find(start_name);
    }
    return start_node;
  }

  int checkEndState() {
    int end_node = -1;
    if (end_name.empty()) {
      vector<size_t> n_out_arcs(nodes.size(), 0);
      for (vector<ThotLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
        n_out_arcs [link->start]++;
      }
      for (size_t n = 0; n < nodes.size(); n++) {
        if (n_out_arcs[n] == 0) {
          if (end_node == -1) end_node = n;
          else {
            LOG(FAIL) << " invalid lattice '" << source << "'; it has more than one end state ";
          }
        }
      }
      if (end_node != -1) end_name = nodes[end_node].name;
    }
    else {
      end_node = ssyms->Find(end_name);
    }
    return end_node;
  }

  template<typename Arc>
  fst::MutableFst<Arc> *CreateFst() {
    fst::VectorFst<Arc> *fst = new fst::VectorFst<Arc>();


    if (checkStartState() == -1) {
      LOG(FAIL) << " invalid lattice '" << source << "'; it does not have initial state ";
    }

    if (checkEndState() == -1) {
      LOG(FAIL) << " invalid lattice '" << source << "'; it does not have final state ";
    }

    for (vector<ThotNode>::const_iterator node = nodes.begin(); node < nodes.end(); ++node) {
      fst->AddState();
    }
    fst->SetStart(checkStartState());
    fst->SetFinal(checkEndState(), Arc::Weight::One());


    for (vector<ThotLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
      Arc arc;
      arc.ilabel = link->words;
      arc.olabel = arc.ilabel;
      arc.weight = compute_weight<typename Arc::Weight>()(*link, weights);
      arc.nextstate = link->end;
      fst->AddArc(link->start, arc);
    }

    fst->SetInputSymbols(wsyms);
    fst->SetOutputSymbols(wsyms);

    return fst;
  }

  template<typename Arc>
  Lattice<Arc> *CreateLattice() {
    Lattice<Arc> *lat = new Lattice<Arc>();

    lat->setWeights(weights.getWeights()); 
    lat->setFeatureNames(weights.getFeatureNames()); 

    if (checkStartState() == -1) {
      LOG(FAIL) << " invalid lattice '" << source << "'; it does not have initial state ";
    }

    if (checkEndState() == -1) {
      LOG(FAIL) << " invalid lattice '" << source << "'; it does not have final state ";
    }

    for (vector<ThotNode>::const_iterator node = nodes.begin(); node < nodes.end(); ++node) {
      lat->getFst().AddState();
    }
    lat->getFst().SetStart(checkStartState());
    lat->getFst().SetFinal(checkEndState(), Arc::Weight::One());

    for (vector<ThotLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
      Arc arc;
      arc.ilabel = link->words;
      arc.olabel = arc.ilabel;
      arc.weight = compute_weight<typename Arc::Weight>()(*link, weights);
      arc.nextstate = link->end;
      lat->getFst().AddArc(link->start, arc);
    }

    lat->getFst().SetInputSymbols(wsyms);
    lat->getFst().SetOutputSymbols(wsyms);

    return lat;
  }

};



class ThotContext {
public:

   void* scanner;   // the scanner state
   std::deque<int> stack; // stack for scanner context nesting

   istream& is;     // input stream
   string source;   // filename

   ThotLattice thot;
   ThotLink *link_ptr;
   ThotNode *node_ptr;

   size_t num_nodes;
   size_t num_links;


   ThotContext(istream& _is = cin, const string &_source = "stdin", const Wordlist& epsilon_symbols = Wordlist()): scanner(0), is(_is), source(_source),
       thot(epsilon_symbols), link_ptr(0), node_ptr(0), num_nodes(0), num_links(0)
   {
//      init_scanner();
      thot.source = source;
   }

   virtual ~ThotContext() {
//      destroy_scanner();
   }

//// Defined in thot.l
//protected:
//   void init_scanner();
//   void destroy_scanner();
};


//int thot_parse(ThotContext*);

} }


//// defines for parser and scanner
//#include "yacc.hpp"
//#include "location.hh"
//#include "position.hh"
//#include "stack.hh"
//#include "../parser-utils.h"
//#include <openlat/iofilter.h>
//
//#define YYSTYPE openlat::thot::parser::semantic_type
//#define YYLTYPE openlat::thot::location
//#define YY_EXTRA_TYPE openlat::thot::ThotContext*
///*  #define YY_USER_ACTION yylloc->first_line = yylineno; */
//
//namespace openlat { namespace thot {
//void thot_error(YYLTYPE* locp, ThotContext* context, const char* err);
//} }
//
//#define YY_INPUT(buf,result,max_size) result = istream_input<ThotContext>(*yyextra, buf, max_size)



#endif /* openlat_THOT_H_ */
