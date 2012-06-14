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
 * htk-printer.h
 *
 *  Created on: 15/03/2012
 *      Author: valabau
 */

#ifndef openlat_HTK_PRINTER_H_
#define openlat_HTK_PRINTER_H_

#include <fst/shortest-distance.h>

namespace openlat {

template<class Arc>
float HtkWeight(const typename Arc::Weight &weight) {
  return -weight.Value();
}

template<class Arc>
void PrintHtk(const fst::Fst<Arc> &fst, std::ostream &out) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  out << "VERSION=1.0\n";

  size_t num_states = 0;
  size_t num_arcs = 0;

  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    num_states++;
    num_arcs += fst.NumArcs(siter.Value());
  }

  out << "N=" << num_states << " L=" << num_arcs << "\n";

  size_t i = 0;
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    out << "I=" << (i++) << "\n";
  }

  const fst::SymbolTable *syms = fst.InputSymbols();

  size_t j = 0;
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();
    for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      const Arc &arc = aiter.Value();
      out << "J=" << (j++) << " S=" << s << " E=" << arc.nextstate << " W=";

      if (syms != 0) {
        out << syms->Find(arc.ilabel);
      }
      else {
        out << arc.ilabel;
      }

      out << " l=" << HtkWeight<Arc>(arc.weight) << "\n";
    }
  }

}

}

#endif /* openlat_HTK_PRINTER_H_ */
