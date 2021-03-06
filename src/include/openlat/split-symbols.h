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
 * split-symbols.h
 *
 *  Created on: 12/04/2012
 *      Author: valabau
 */

#ifndef openlat_SPLIT_SYMBOLS_H_
#define openlat_SPLIT_SYMBOLS_H_

#include <fst/mutable-fst.h>
#include <unordered_map>
#include <boost/shared_ptr.hpp>

namespace openlat {

template <typename A> struct TrieNode;
template <typename A> struct TrieNode: std::map< typename A::Label, std::pair< typename A::StateId, boost::shared_ptr< TrieNode<A> > > > {};


template<class Arc>
void SplitSymbols(const fst::Fst<Arc> &fst, fst::MutableFst<Arc> *ofst,
                  const std::string &separator = "", const std::string *space = 0,
                  const fst::SymbolTable *const_syms = 0)
{
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::Label Label;
  typedef typename TrieNode<Arc>::value_type TrieNodeValue;
  typedef typename TrieNode<Arc>::iterator   TrieNodeIter;

  bool is_acceptor = fst.Properties(fst::kAcceptor, true);
//  bool is_deterministic = fst.Properties(fst::kIDeterministic, true);
  fst::SymbolTable isyms("inputs");
  fst::SymbolTable osyms("outputs");

  isyms.AddSymbol(fst.InputSymbols()->Find(int64(0)));
  osyms.AddSymbol(fst.OutputSymbols()->Find(int64(0)));


  typedef std::unordered_map<StateId, StateId> StateMap;
  StateMap state_map;

  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId state = ofst->AddState();
    state_map[siter.Value()] = state;


    if (siter.Value() == fst.Start()) ofst->SetStart(state);
  }

  StateId superfinal = ofst->AddState();
  ofst->SetFinal(superfinal, Weight::One());



  std::vector<std::string> tokens;
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    TrieNode<Arc> trie;
    // initialize trie with no symbol to go to start state
    TrieNodeIter root = trie.insert(std::make_pair(fst::SymbolTable::kNoSymbol, std::make_pair(state_map[siter.Value()], boost::shared_ptr< TrieNode<Arc> >(new TrieNode<Arc>)))).first;

    // iterate for each arc in the node
    for (fst::ArcIterator < fst::Fst<Arc> > aiter(fst, siter.Value()); !aiter.Done(); aiter.Next()) {
      const Arc &arc = aiter.Value();
      tokens.clear();

      const std::string &istr = fst.InputSymbols()->Find(arc.ilabel);
      const std::string &ostr = fst.OutputSymbols()->Find(arc.olabel);

      Weight final_weight = fst.Final(arc.nextstate);
      bool is_final = (final_weight.Member() and final_weight != Weight::Zero());

      // the symbol is a constant symbols and should not be split
      if (const_syms != 0 and const_syms->Find(istr) != fst::SymbolTable::kNoSymbol) {
        tokens.push_back(istr);
      }
      // split the string into tokens
      else {
        tokenize(istr, tokens, separator);
        // if tokens is empty just add the original string
        if (tokens.empty()) tokens.push_back(istr);
      }

      // add nodes to the fst
      {
        TrieNodeIter start = root;

        // iterate all tokens but the last one
        for (size_t i = 0; i < tokens.size() - 1; i++) {
          Label label = isyms.AddSymbol(tokens[i]);
          StateId start_state = start->second.first;

          boost::shared_ptr< TrieNode<Arc> > node = start->second.second;
          TrieNodeIter end = node->find(label);
          if (end == node->end()) {
            StateId end_state = ofst->AddState();
            end = node->insert(std::make_pair(label, std::make_pair(end_state, boost::shared_ptr< TrieNode<Arc> >(new TrieNode<Arc>)))).first;

            Arc oarc;
            oarc.ilabel = label;
            oarc.olabel = (is_acceptor)?osyms.AddSymbol(tokens[i]):fst::SymbolTable::kNoSymbol;
            oarc.weight = Weight::One();
            oarc.nextstate = end_state;
            ofst->AddArc(start_state, oarc);
          }

          start = end;
        }

        // add last token

        // if the original state was final add arc to the superfinal state
        if (is_final) {
          StateId start_state = start->second.first;

          Arc oarc;
          oarc.ilabel = isyms.AddSymbol(tokens.back());
          oarc.olabel = (is_acceptor)?osyms.AddSymbol(tokens.back()):osyms.AddSymbol(ostr);
          oarc.weight = Times(arc.weight, final_weight);
          oarc.nextstate = superfinal;
          ofst->AddArc(start_state, oarc);
        }

        // if next state has outgoing arcs
        fst::ArcIterator < fst::Fst<Arc> > nsaiter(fst, arc.nextstate);
        if (not nsaiter.Done()) {
          // add last token
          StateId start_state = start->second.first;
          StateId end_state = state_map[arc.nextstate];
          if (space != 0) {
            end_state = ofst->AddState();
          }

          Arc oarc;
          oarc.ilabel = isyms.AddSymbol(tokens.back());
          oarc.olabel = (is_acceptor)?osyms.AddSymbol(tokens.back()):osyms.AddSymbol(ostr);
          oarc.weight = arc.weight;
          oarc.nextstate = end_state;
          ofst->AddArc(start_state, oarc);

          // add space token
          if (space != 0) {
            start_state = end_state;
            end_state = state_map[arc.nextstate];

            oarc.ilabel = isyms.AddSymbol(*space);
            oarc.olabel = (is_acceptor)?osyms.AddSymbol(*space):fst::SymbolTable::kNoSymbol;
            oarc.weight = Weight::One();
            oarc.nextstate = end_state;
            ofst->AddArc(start_state, oarc);
          }
        }

      }
    }
  }
  Connect(ofst);

  ofst->SetInputSymbols(&isyms);
  ofst->SetOutputSymbols(&osyms);
}

}

#endif /* openlat_SPLIT_SYMBOLS_H_ */
