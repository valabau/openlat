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
 * thot-printer.h
 *
 *  Created on: 15/03/2012
 *      Author: valabau
 */

#ifndef openlat_THOT_PRINTER_H_
#define openlat_THOT_PRINTER_H_

#include <fst/shortest-distance.h>

namespace openlat {

template<class Arc>
float ThotWeight(const typename Arc::Weight &weight) {
  return -weight.Value();
}

template<class Arc>
typename Arc::StateId SwapStateId(const typename Arc::StateId id, const typename Arc::StateId id1, const typename Arc::StateId id2) {
  if (id == id1) return id2;
  else if (id == id2) return id1;
  return id;
}


template<class Arc>
void PrintThot(const fst::Fst<Arc> &fst, std::ostream &out) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  out << "# Generated by openlat thot-printer\n";

  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    if (fst.Final(siter.Value()) != Weight::Zero()) { out << siter.Value() << " "; }
  }
  out << "\n";

  const fst::SymbolTable *syms = fst.InputSymbols();

  // print first start state and replace state 0 with start id
  StateId start = fst.Start();
  StateId idzero = 0;
  {
    for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, start); !aiter.Done(); aiter.Next()) {
      const Arc &arc = aiter.Value();
      out << idzero << " " << SwapStateId<Arc>(arc.nextstate, idzero, start) << " " << ThotWeight<Arc>(arc.weight) << " ";

      if (arc.ilabel == 0) {
        out << "!NULL";
      }
      else if (syms != 0) {
        out << syms->Find(arc.ilabel);
      }
      else {
        out << arc.ilabel;
      }
      out << "\n";
    }
  }


  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    if (siter.Value() != start) {
      StateId s = SwapStateId<Arc>(siter.Value(), idzero, start);

      for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, siter.Value()); !aiter.Done(); aiter.Next()) {
        const Arc &arc = aiter.Value();
        out << s << " " << SwapStateId<Arc>(arc.nextstate, idzero, start) << " " << ThotWeight<Arc>(arc.weight) << " ";

        if (arc.ilabel == 0) {
          out << "!NULL";
        }
        else if (syms != 0) {
          out << syms->Find(arc.ilabel);
        }
        else {
          out << arc.ilabel;
        }
        out << "\n";
      }
    }
  }

}

template<class Arc>
void AddStartAndEndSymbols(fst::MutableFst<Arc> *fst, const char *start_symbol = 0, const char *end_symbol = 0) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  if (start_symbol == 0 and end_symbol == 0) return;

  fst::SymbolTable *isyms = fst->MutableInputSymbols();
  fst::SymbolTable *osyms = fst->MutableOutputSymbols();

  if (start_symbol) {
    typename Arc::Label istart = isyms->AddSymbol(start_symbol);
    typename Arc::Label ostart = osyms->AddSymbol(start_symbol);

    StateId old_start = fst->Start();
    StateId new_start = fst->AddState();
    fst->SetStart(new_start);
    fst->AddArc(new_start, Arc(istart, ostart, Weight::One(), old_start));
  }

  if (end_symbol) {
    typename Arc::Label iend = isyms->AddSymbol(end_symbol);
    typename Arc::Label oend = osyms->AddSymbol(end_symbol);

    StateId new_end = fst->AddState();
    fst->SetFinal(new_end, Weight::One());

    for (fst::StateIterator<fst::Fst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
      StateId s = siter.Value();
      if (fst->Final(s) != Weight::Zero()) {
        fst->AddArc(s, Arc(iend, oend, fst->Final(s), new_end));
        fst->SetFinal(s, Weight::Zero());
      }
    }
  }
}

}

#endif /* openlat_THOT_PRINTER_H_ */
