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
 * htk.h
 *
 *  Created on: 23/02/2012
 *      Author: valabau
 */

#ifndef openlat_HTK_H_
#define openlat_HTK_H_

#include <string>
#include <vector>
#include <deque>
#include <tr1/unordered_map>

#include <fst/fst.h>
#include <fst/arc.h>
#include <fst/symbol-table.h>
#include <openlat/utils.h>
#include <openlat/vector-weight.h>
#include <openlat/htk-compiler.h>
#include <cmath>

using namespace std;
using namespace fst;


namespace openlat { namespace htk {

class HtkLink {
public:
  string name;     // Link name
  int start; // Start node number (of the link)
  int end;   // End node number (of the link)
  int input;    // Input symbol
  int output;   // Output symbol
  int var;         // Pronunciation variant number
  string div;      // Segmentation [ d=:(label[,duration,like]:)+ ]
  vector<float> features; // Features

  HtkLink(): start(kNoStateId), end(kNoStateId), input(kNoLabel), output(kNoLabel), var(0) {};
};

class HtkNode {
public:
  string name;     // Node name
  float time;      // Time from start of utterance (in seconds)
  string sublat;   // Substitute named sub-lattice for this node

  HtkLink link;    // auxiliar link
  HtkNode() {};
};

class HtkWeights {
public:
  typedef enum {ACOUSTIC, LANGUAGE, NGRAM, LMIN, LMOUT, WDPENALTY, OUTPUT_WDPENALTY, NOISE_PENALTY, POSTERIOR, XSCORE} feature_t;

  struct feature_hash {
    size_t operator()(const feature_t& __x) const {
      return std::tr1::hash<int>()(static_cast<int>(__x));
    }
  };

  HtkWeights(): amscale(1.0) {}

  vector<float> weights; // scales for features
  float  amscale;  // Acoustic scale for posterior computation
  unordered_map<feature_t, size_t, feature_hash> feature_pos; // position of the feature

  size_t getFeature(feature_t feat_name) {
    unordered_map<feature_t, size_t>::const_iterator it = feature_pos.find(feat_name);
    if (it == feature_pos.end()) {
      it = feature_pos.insert(make_pair(feat_name, feature_pos.size())).first;
    }
    return it->second;
  }

  std::vector<std::string> getFeatureNames() {
    unordered_map<feature_t, size_t>::const_iterator it;
    std::string strnames[] = {"acoustic", "language", "ngram", "lmin", "lmout", "wdpenalty", "output_wdpenalty", "noise_penalty", "posterior", "xscore" };

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
  inline Weight operator()(const HtkLink &link, const HtkWeights &weights) {
    float w = 0;
    const size_t m = std::min(link.features.size(), weights.weights.size());
    size_t i;
    for (i = 0; i < m; i++) {
      w += link.features[i] * weights.weights[i];
    }
    for (; i < link.features.size(); i++) {
      w += link.features[i];
    }
    w /= -weights.amscale;
    return Weight(w);
  }
};

template <>
struct compute_weight<LogLinearWeight> {
  inline LogLinearWeight operator()(const HtkLink &link, const HtkWeights &weights) {
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



class HtkLattice {
public:

  string version;   // Lattice specification adhered to
  string utterance; // Utterance identifier
  string sublat;    // Sub-lattice name
  float  base;      // LogBase for Likelihoods (0.0 not logs,default base e)
  string vocab;     // name of the vocabulary
  string lmname;    // Name of Language model
  string source;    // filename


  SymbolTable *ssyms;
  SymbolTable *asyms;
  SymbolTable *isyms;
  SymbolTable *osyms;

  HtkWeights weights;
  vector<HtkLink> links;
  vector<HtkNode> nodes;

  string start_name;
  string end_name;

  bool has_output;

  HtkLattice(): base(M_E),
                ssyms(new SymbolTable("states")), asyms(new SymbolTable("arcs")),
                isyms(new SymbolTable("inputs")), osyms(new SymbolTable("outputs")),
                has_output(false)
  {
    // add epsilon symbols (!NULL in htk)
    isyms->AddSymbol("!NULL");
    osyms->AddSymbol("!NULL");
    isyms->AddSymbol("<s>");
    osyms->AddSymbol("<s>");
    isyms->AddSymbol("</s>");
    osyms->AddSymbol("</s>");
  };

  ~HtkLattice() {
    delete ssyms; delete asyms; delete isyms; delete osyms;
  };

  float getScore(const vector<float> scores) {
    float sum = .0;
    if (scores.size() != weights.weights.size()) LOG(FATAL) << "Score size and weigth size do not match\n";
    for (size_t i = 0; i < scores.size(); i++) {
      sum += scores[i] * weights.weights[i];
    }
    return sum / weights.amscale;
  }

  HtkLink *getArc(const string arc_name) {
    size_t id = asyms->AddSymbol(arc_name);
    if (links.size() < id + 1) links.resize(id + 1);
    return &links[id];
  }

  HtkNode *getState(string state_name) {
    size_t id = ssyms->AddSymbol(state_name);
    if (nodes.size() < id + 1) nodes.resize(id + 1);
    nodes[id].name = state_name;
    nodes[id].link.name = state_name + "_" + "aux";
    return &nodes[id];
  }

  LogLinearFst::StateId getStateId(string state_name) {
    return ssyms->AddSymbol(state_name);
  }



  void assign(HtkLink *l, HtkWeights::feature_t feat_name, float value) {
    size_t pos = weights.getFeature(feat_name);
    if (pos >= l->features.size())     l->features.resize(pos + 1, .0);
    if (pos >= weights.weights.size()) weights.weights.resize(pos + 1, 1.0);
    // change log base to base e and negate for the log semiring
    l->features[pos] = change_base(value, base);
  }

  void assignWeight(HtkWeights::feature_t feat_name, float value) {
    size_t pos = weights.getFeature(feat_name);
    if (pos >= weights.weights.size()) weights.weights.resize(pos + 1, 1.0);
    weights.weights[pos] = value;
  }

  void assignInput(HtkLink *l, const string& input) {
    l->input = isyms->AddSymbol(input);
  }

  void assignOutput(HtkLink *l, const string& output) {
    has_output = true;
    l->output = osyms->AddSymbol(output);
  }

  int checkStartState() {
    int start_node = -1;
    if (start_name.empty()) {
      vector<size_t> n_in_arcs(nodes.size(), 0);
      for (vector<HtkLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
        n_in_arcs [link->end]++;
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
      for (vector<HtkLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
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

    for (vector<HtkNode>::const_iterator node = nodes.begin(); node < nodes.end(); ++node) {
      fst->AddState();
    }
    fst->SetStart(checkStartState());
    fst->SetFinal(checkEndState(), Arc::Weight::One());


    for (vector<HtkLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
      Arc arc;
      arc.ilabel = link->input;
      arc.olabel = ((has_output)?link->output:link->input);
      arc.weight = compute_weight<typename Arc::Weight>()(*link, weights);
      arc.nextstate = link->end;
      fst->AddArc(link->start, arc);
    }

    fst->SetInputSymbols(isyms);
    if (has_output) {
      fst->SetOutputSymbols(osyms);
    }
    else {
      fst->SetOutputSymbols(isyms);
    }

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

    for (vector<HtkNode>::const_iterator node = nodes.begin(); node < nodes.end(); ++node) {
      lat->getFst().AddState();
    }
    lat->getFst().SetStart(checkStartState());
    lat->getFst().SetFinal(checkEndState(), Arc::Weight::One());

    for (vector<HtkLink>::const_iterator link = links.begin(); link < links.end(); ++link) {
      Arc arc;
      arc.ilabel = link->input;
      arc.olabel = ((has_output)?link->output:link->input);
      arc.weight = compute_weight<typename Arc::Weight>()(*link, weights);
      arc.nextstate = link->end;
      lat->getFst().AddArc(link->start, arc);
    }

    lat->getFst().SetInputSymbols(isyms);
    if (has_output) {
      lat->getFst().SetOutputSymbols(osyms);
    }
    else {
      lat->getFst().SetOutputSymbols(isyms);
    }

    return lat;
  }

};



class HtkContext {
public:

   void* scanner;   // the scanner state
   std::deque<int> stack; // stack for scanner context nesting

   istream& is;     // input stream
   string source;   // filename

   HtkLattice htk;
   HtkLink *link_ptr;
   HtkNode *node_ptr;

   size_t num_nodes;
   size_t num_links;


   HtkContext(istream& _is = cin, const string &_source = "stdin"): scanner(0), is(_is), source(_source),
       link_ptr(0), node_ptr(0), num_nodes(0), num_links(0)
   {
      init_scanner();
      htk.source = source;
   }

   virtual ~HtkContext() {
      destroy_scanner();
   }

// Defined in htk.l
protected:
   void init_scanner();
   void destroy_scanner();
};


int htk_parse(HtkContext*);

} }


// defines for parser and scanner
#include "yacc.hpp"
#include "location.hh"
#include "position.hh"
#include "stack.hh"
#include "../parser-utils.h"
#include <openlat/iofilter.h>

#define YYSTYPE openlat::htk::parser::semantic_type
#define YYLTYPE openlat::htk::location
#define YY_EXTRA_TYPE openlat::htk::HtkContext*
/*  #define YY_USER_ACTION yylloc->first_line = yylineno; */

namespace openlat { namespace htk {
void htk_error(YYLTYPE* locp, HtkContext* context, const char* err);
} }

#define YY_INPUT(buf,result,max_size) result = istream_input<HtkContext>(*yyextra, buf, max_size)



#endif /* openlat_HTK_H_ */
