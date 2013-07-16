/*
 *   Copyright 2013, valabau
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
 * confusion-network.h
 *
 *  Created on: 09/05/2013
 *      Author: valabau
 */

#ifndef OPENLAT_CONFUSION_NETWORK_H_
#define OPENLAT_CONFUSION_NETWORK_H_

#include <fst/shortest-distance.h>
#include <openlat/normalize.h>

namespace openlat {

//// Confision Network auxiliar structures
//template<class Arc>
//class CnNode {
//public:
//  typedef std::map<typename Arc::Label, std::vector<const Arc*> > words_t;
//  typedef std::set<const Arc*> arcs_t;
//  float time;
//  words_t words;
//  arcs_t arcs;
//  CnNode(float _time): time(_time)  { };
//  void addLink(const Arc* link, const typename Arc::Label &word) {
//    arcs.insert(link);
//    typename words_t::const_iterator it = words.find(word);
//    if (it == words.end()) it = words.insert(make_pair(word, vector<const Arc*>())).first;
//    it->second.push_back(link);
//  };
//  bool hasPredecessorOf(size_t) { return false; }
//};
//
//
///*
// * Creates a list of cn_nodes_t representing a confusion network
// * using the pivot algorithm (D. Hakkani-Tür et al. / Computer Speech and Language 20 (2006) 495-514)
// * Link features are still to be merged
// */
//template<class Arc>
//void ComputeConfusionNodes(const fst::Fst<Arc> &fst, std::list< CnNode<Arc> >& cn_nodes) {
//  // 1. Compute word posterior probabilities
//  typedef typename Arc::Weight Weight;
//  typedef typename Arc::StateId StateId;
//
//  std::vector<Weight> forward, backward;
//  fst::ShortestDistance(fst, &forward, false);
//  fst::ShortestDistance(fst, &backward, true);
//
//  // 2. Extract the pivot baseline path
//  {
//    std::vector<const Link *> path;
//    vector<W> probs;
//    // perform a standard viterbi decoding using recognition scores
//    decodeViterbi(path, probs, DECODING_MODE_RECOGNITION, initial);
//
//    assert_bt(not path.empty(), "Decode Viterbi could not found a complete hypothesis");
//    for (size_t i = 0; i < path.size(); i++) {
//      CnNode node(nodes[path[i]->start]->time);
//      cn_nodes.push_back(node);
//    }
//    {
//      CnNode node(nodes[path.back()->end]->time);
//      cn_nodes.push_back(node);
//    }
//  }
//
//  // Compute predecessor matrix.
//  vector< set<link_id_t> > predecessors(nodes.size());
//  for (vector<Node *>::const_iterator nit = nodes.begin(); nit < nodes.end(); ++nit) {
//    for (vector<Link *>::const_iterator lit = (*nit)->incoming.begin();
//         lit < (*nit)->incoming.end(); ++lit)
//    { // store just first order predecessors
//      predecessors[(*nit)->id].insert((*lit)->id);
//    }
//  }
//  // find the last node, as it will remain the last node during the whole the process
//  list<CnNode>::const_iterator final_node_it = cn_nodes.end();
//  final_node_it--;
//
//   // 3. For all links l in topological order
//  // For the algorithm to work correctly with only first order predecessors
//  // we must iterate over the links sorted by the central time they occur
//  for (size_t l = 0; l < arcs.size(); l++) {
//    const Link * link = arcs[l];
//    // 3.1. Find the two consecutive states S_s and S_e that have the maximum overlap,
//    //      O(l, S_s, S_e), with l
//    W max_overlap = - numeric_limits<W>::infinity();
//    list<CnNode>::iterator max_it = cn_nodes.end();
//
//    // we take advantadge of the fact that links are sorted so the maximum overlap
//    // happens either with current node or with the next node
//    for (list<CnNode>::iterator it = cn_nodes.begin(); it != final_node_it; ++it)
//    {
//      list<CnNode>::iterator nit = it; nit++;
//      W overlap = min(nit->time, link->endTime()) - max(it->time, link->startTime());
//      if (overlap > max_overlap) {
//        max_overlap = overlap;
//        max_it = it;
//      }
//    }
//
//    // 3.2. If there is no transition in between S_s and S_e that precedes l in a path
//    //      of the lattice, insert l between them
//    set<link_id_t> s;
//    set_intersection(predecessors[link->start].begin(), predecessors[link->start].end(),
//                     max_it->arcs.begin(), max_it->arcs.end(),
//                     inserter(s, s.begin()));
//
//    WordId word_id = arc.ilabel;
//    if (not wordmap.getWord(arc.ilabel).isEvent()) word_id = wordmap.getWordId(w_eps);
//
//    if (s.empty()) {
//      max_it->addLink(static_cast<link_id_t>(l), word_id);
//    }
//    // 3.3. Otherwise,
//    else {
//      // 3.3.a. Create a new state S_n
//      float ts = max_it->time;
//      max_it++;
//      float te = max_it->time;
//      // we assign the time as the mean of ts and te
//      struct CnNode node((ts + te) / 2);
//      max_it = cn_nodes.insert(max_it, node);
//      // 3.3.b. Change the destination state of all transitions originating from state S_s to S_n
//      //        * Not necessary in this implementation
//      // 3.3.c. Insert l between S_n and S_e
//      max_it->addLink(static_cast<link_id_t>(l), word_id);
//    }
//
//  }
//
//
//  // 4. Insert epsilon links in between every two consecutive states
//  // This must be done at a latter stage
//
//}
//
///*
// * Creates a Confusion Network
// * Creates a Confusion Network using the pivot algorithm (D. Hakkani-Tür et al. / Computer Speech and Language 20 (2006) 495-514)
// */
//template<class Arc>
//void ConfusionNetwork(const fst::Fst<Arc> &fst, fst::MutableFst<Arc> *ofst)
//{
//  typedef typename Arc::Weight Weight;
//  typedef typename Arc::StateId StateId;
//
//  std::vector<Weight> forward, backward;
//  fst::ShortestDistance(fst, &forward, false);
//  fst::ShortestDistance(fst, &backward, true);
//
//  vector<CnNode> cn_nodes;
//  {
//    list<CnNode> list_cn_nodes;
//    ComputeConfusionNodes(fst, list_cn_nodes);
//
//    size_t n = 0, n_id = 0;
//    size_t last_cn_node = list_cn_nodes.size() - 1;
//    for (list<CnNode>::iterator nit = list_cn_nodes.begin(); n < last_cn_node; ++nit, n++) {
//      Weight sum = Weight::Zero();
//      for (CnNode::words_t::const_iterator wit = nit->words.begin();
//           wit != nit->words.end(); ++wit)
//      {
//        if (wit->first != 0) { // the word is not eps
//          Weight word_sum = Weight::Zero();
//          for (std::vector<const Arc*>::const_iterator slit = wit->second.begin();
//               slit != wit->second.end(); ++slit)
//          {
//            const Arc& arc = *(*slit);
//            Weight posterior = fst::Plus(forward[arc.nextstate], fst::Plus(arc.weight, backward[arc.nextstate]));
//            word_sum = fst::Plus(word_sum, posterior);
//          }
//          sum = fst::Plus(sum, word_sum);
//        }
//      }
//
//      // if posterior is 0 we can remove the node
//      // else add node
//      if (sum != Weight::Zero()) {
//        cn_nodes.push_back(*nit);
//        ofst->AddState();
//      }
//    }
//
//    // add final node
//    cn_nodes.push_back(list_cn_nodes.back());
//    StateId final = ofst->AddState();
//    ofst->SetFinal(final, Weight::One());
//  }
//
//  // copy to new slf
//  {
//    size_t n = 0;
//    for (vector<CnNode>::const_iterator nit = cn_nodes.begin(); nit < cn_nodes.end() - 1; ++nit, n++)
//    {
//      if (not nit->arcs.empty()) {
//        node_id_t sn = slf.getInternalNodeId(n);
//        node_id_t en = slf.nodes[sn + 1]->id;
//        //node_id_t en = slf.getInternalNodeId(n + 1);
//        assert_bt(sn != UNK_NODE and en != UNK_NODE, "Invalid CN topology");
//
//        for (CnNode::words_t::const_iterator wit = nit->words.begin();
//             wit != nit->words.end(); ++wit)
//        {
//          W posterior = LOG_ZERO;
//          for (vector<link_id_t>::const_iterator slit = wit->second.begin();
//               slit != wit->second.end(); ++slit)
//          {
//            posterior = add_log(posterior, link_post[*slit]);
//          }
//
//          if (posterior != Weight::Zero()) {
//            Link * link_ptr     = new Link(slf.features.size());
//            link_ptr->start     = sn;
//            link_ptr->end       = en;
//            link_ptr->word_id   = slf.wordmap.getWordId(wit->first.getWord());
//            slf.addLink(link_ptr);
//            link_ptr->posterior = posterior;
//
//            eps_post = add_log(eps_post, posterior);
//          }
//        }
//
//        // Insert epsilon
//        if (not nit->words.empty()) {
//          bool has_eps = (eps_link != NULL);
//          if (not has_eps) {
//            Link * link_ptr = new Link(slf.features.size());
//            link_ptr->start = sn;
//            link_ptr->end   = en;
//            link_ptr->word_id   = slf.wordmap.getWordId(w_eps);
//            slf.addLink(link_ptr);
//            eps_link = link_ptr;
//          }
//          eps_link->posterior = sub_log(W(LOG_ONE), eps_post);
//        }
//      }
//    delete[] eps_feats;
//  }
//  delete[] link_post;
//}


template <typename Arc>
void ComputeMinimumBayesRiskTable(const fst::Fst<Arc> &fst, std::vector< std::map<typename Arc::Label, typename Arc::Weight> > &table) {
  typedef typename Arc::Weight W;
  typedef typename Arc::StateId StateId;

  std::vector<W> backward;
  fst::ShortestDistance(fst, &backward, true);

  typedef std::map<size_t, W> PosMap;
  typedef std::vector<PosMap> StatePosMap;
  StatePosMap forward(NumStates(fst));
  forward[fst.Start()][0] = W::One(); 


  /* for each node in order  */
  for (fst::StateIterator<fst::Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    StateId s = siter.Value();

    /* for each incoming link */
    for (fst::ArcIterator< fst::Fst<Arc> > aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();

      /* iterate over all lengths */
      for (typename PosMap::const_iterator mit = forward[s].begin(); mit != forward[s].end(); ++mit) {

        // get position of current word and resize table if necessary
        size_t position = mit->first;
        if (arc.ilabel != 0) {
          position++;
        }
        if (table.size() < position + 1) table.resize(position + 1);

        // add up \Phi^(i) forward probabilities for position (i)
        W fwd = fst::Times(mit->second, arc.weight);
        if (forward[arc.nextstate].find(position) != forward[arc.nextstate].end()) {
          forward[arc.nextstate][position] = fst::Plus(forward[arc.nextstate][position], fwd);
        }
        else {
          forward[arc.nextstate].insert(make_pair(position, fwd));
        }

        // add up word posterior probabilities for word in position (i)
        if (arc.ilabel != 0) {
          W post = fst::Times(fwd, backward[arc.nextstate]);
          if (table[position].find(arc.ilabel) != table[position].end()) {
            table[position][arc.ilabel] = fst::Plus(table[position][arc.ilabel], post);
          } else {
            table[position].insert(make_pair(arc.ilabel, post));
          }
        }
      }
    }

    if (fst.Final(s) != W::Zero()) {
      for (typename PosMap::const_iterator mit = forward[s].begin(); mit != forward[s].end(); ++mit) {
        size_t position = mit->first + 1;
        if (table.size() < position + 1) table.resize(position + 1);
        W post = fst::Times(mit->second, backward[s]);
        typename Arc::Label ilabel = fst::kNoLabel;
        if (table[position].find(ilabel) != table[position].end()) {
          table[position][ilabel] = fst::Plus(table[position][ilabel], post);
        } else {
          table[position].insert(make_pair(ilabel, post));
        }
      }
    }

  }

  typedef std::map<typename Arc::Label, W> WordMap;
  typedef std::vector<WordMap> Table;
  size_t s = 0;
  for (typename Table::iterator tit = table.begin(); tit < table.end(); ++tit, ++s) {
    W sum = W::Zero();
    for (typename WordMap::const_iterator it = tit->begin(); it != tit->end(); ++it) {
      sum = fst::Plus(sum, it->second);
    }
    for (typename WordMap::iterator it = tit->begin(); it != tit->end(); ++it) {
      it->second = fst::Divide(it->second, sum);
    }
  }

}



/*
 * Creates a Position Word Posterior lattice
 *
 */
template <typename Arc>
void PositionPosteriorNetwork(const fst::Fst<Arc> &fst, fst::MutableFst<Arc> *ofst) {
  typedef typename Arc::Weight W;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;

  typedef std::map<typename Arc::Label, W> WordMap;
  typedef std::vector<WordMap> Table;

  Table table;
  ComputeMinimumBayesRiskTable(fst, table);

  ofst->DeleteStates();
  ofst->SetInputSymbols(fst.InputSymbols());
  ofst->SetOutputSymbols(fst.OutputSymbols());

  ofst->SetStart(0);
  for (size_t s = 0; s < table.size(); s++) ofst->AddState();
  ofst->SetFinal(table.size() - 1, W::One());

  for (size_t s = 0; s < table.size(); s++) {
    if (not table[s].empty()) {
      for (typename WordMap::const_iterator it = table[s].begin(); it != table[s].end(); ++it) {
        //std::string  str = fst.InputSymbols()->Find(it->first);
        //std::cerr << s << ": " << str << " = " << exp(-it->second.Value()) << std::endl;

        if (it->first == fst::kNoLabel) {
          ofst->SetFinal(s - 1, it->second);
        }
        else {
          ofst->AddArc(s - 1, Arc(it->first, it->first, it->second, s));
        }
      }
    }
  }

  fst::Connect(ofst);
}



}

#endif /* OPENLAT_CONFUSION_NETWORK_H_ */
