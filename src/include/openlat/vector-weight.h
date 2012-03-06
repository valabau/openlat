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
 * vector-weight.h
 *
 *  Created on: 21/02/2012
 *      Author: valabau
 */

#ifndef openlat_VECTOR_WEIGHT_H_
#define openlat_VECTOR_WEIGHT_H_


#include <string>
#include <vector>

#include <fst/weight.h>
#include <fst/vector-fst.h>


DECLARE_string(fst_weight_parentheses);
DECLARE_string(fst_weight_separator);

namespace openlat {

// vector weight
template <class W>
class VectorWeight {
  typedef enum { VECTOR_WEIGHT_ZERO, VECTOR_WEIGHT_ONE, VECTOR_WEIGHT_OTHER } vector_weight_t;

  template<class U> friend bool operator== (const VectorWeight<U> &_w1, const VectorWeight<U> &_w2);
  template<class U> friend bool operator!= (const VectorWeight<U> &_w1, const VectorWeight<U> &_w2);
  template<class U> friend bool ApproxEqual(const VectorWeight<U> &_w1, const VectorWeight<U> &_w2, float delta);

public:
  typedef W Weight;
  typedef VectorWeight<typename W::ReverseWeight> ReverseWeight;

  //  VectorWeight() {}

  VectorWeight(vector_weight_t type = VECTOR_WEIGHT_ONE): type_(type) {}

  VectorWeight(const VectorWeight &w): type_(VECTOR_WEIGHT_OTHER) {
    values_ = w.values_;
  }

  VectorWeight(const vector<W> &w): type_(VECTOR_WEIGHT_OTHER) {
    values_ = w;
  }

  template <class Iterator>
  VectorWeight(Iterator begin, Iterator end): type_(VECTOR_WEIGHT_OTHER) {
    values_.resize(end - begin);
    copy(begin, end, values_.begin());
  }

  VectorWeight(const W &w, size_t n): type_(VECTOR_WEIGHT_OTHER) {
    values_ = vector<W>(n, w);
  }

  void Resize(size_t n) {
    if (type_ == VECTOR_WEIGHT_ZERO) {
      type_ = VECTOR_WEIGHT_OTHER;
      values_.clear();
      values_.resize(n, W::Zero());
    }
    else if (type_ == VECTOR_WEIGHT_ONE) {
      type_ = VECTOR_WEIGHT_OTHER;
      values_.clear();
      values_.resize(n, W::One());
    }
    else {
      values_.resize(n, W::One());
    }
  }

  size_t Length() const {
    return values_.size();
  }


  static const VectorWeight<W> &Zero() {
    static const VectorWeight<W> zero(VECTOR_WEIGHT_ZERO);
    return zero;
  }

  static const VectorWeight<W> &One() {
    static const VectorWeight<W> one(VECTOR_WEIGHT_ONE);
    return one;
  }

  static const ReverseWeight &ReverseZero() {
    static const ReverseWeight zero(VECTOR_WEIGHT_ZERO);
    return zero;
  }

  static const ReverseWeight &ReverseOne() {
    static const ReverseWeight one(VECTOR_WEIGHT_ONE);
    return one;
  }


  std::ostream &Write(std::ostream &strm) const {
    for (size_t i = 0; i < Length(); ++i)
      values_[i].Write(strm);
    return strm;
  }

  VectorWeight<W> &operator=(const VectorWeight<W> &w) {
    type_ = w.type_;
    values_ = w.values_;
    return *this;
  }

  bool Member() const {
    // Zero and One are already memebers
    bool member = true;
    for (size_t i = 0; i < Length(); ++i)
      member = member && values_[i].Member();
    return member;
  }

  size_t Hash() const {
    if (type_ == VECTOR_WEIGHT_ZERO) return size_t(W::Zero().Hash());
    else if (type_ == VECTOR_WEIGHT_ONE) return size_t(W::One().Hash());
    else {
      uint64 hash = 0;
      for (size_t i = 0; i < Length(); ++i)
        hash = 5 * hash + values_[i].Hash();
      return size_t(hash);
    }
  }

  VectorWeight<W> Quantize(float delta = fst::kDelta) const {
    if (type_ == VECTOR_WEIGHT_OTHER) {
      VectorWeight<W> w;
      w.type_ = type_;
      w.values_.resize(Length());
      for (size_t i = 0; i < Length(); ++i)
        w.values_[i] = values_[i].Quantize(delta);
      return w;
    }
    else return *this;
  }

  ReverseWeight Reverse() const {
    if (type_ == VECTOR_WEIGHT_ZERO) {
      return ReverseZero();
    }
    if (type_ == VECTOR_WEIGHT_ONE) {
      return ReverseOne();
    }
    else {
      VectorWeight<W> w;
      w.values_.resize(Length());
      for (size_t i = 0; i < Length(); ++i)
        w.values_[i] = values_[i].Reverse();
      return w;
    }
  }

  const W& Value(size_t i) const {
    if (type_ == VECTOR_WEIGHT_ZERO)     return W::Zero();
    else if (type_ == VECTOR_WEIGHT_ONE) return W::One();
    else if (i < Length()) {
      return values_[i];
    }
    else {
      LOG(ERROR) << " accessing to vector weight outside boundaries ";
      return W::One();
    }
  }

  static const string &Type(void) {
    static const string type = "vector_" + W::Type();
    return type;
  }


 protected:

 private:
  vector<W> values_;
  vector_weight_t type_;

};

#define __VECTOR_WEIGHT_ASSIGN(x,l) ((x.type_ == VectorWeight<W>::VECTOR_WEIGHT_ZERO)?(VectorWeight<W>(W::Zero(), l)):((x.type_ == VectorWeight<W>::VECTOR_WEIGHT_ONE)?(VectorWeight<W>(W::One(), l)):x))

template <class W>
inline bool operator==(const VectorWeight<W> &_w1,
                       const VectorWeight<W> &_w2) {
  if (_w1.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER and _w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) return _w1.type_ == _w2.type_;
  const VectorWeight<W> &w1 = __VECTOR_WEIGHT_ASSIGN(_w1, _w2.Length());
  const VectorWeight<W> &w2 = __VECTOR_WEIGHT_ASSIGN(_w2, _w1.Length());

  if (w1.Length() != w2.Length()) {
    LOG(ERROR) << " comparing vectors of different length ";
    return false;
  }

  bool equal = true;
  for (size_t i = 0; i < w1.Length(); ++i)
    equal = equal && (w1.Value(i) == w2.Value(i));
  return equal;
}

template <class W>
inline bool operator!=(const VectorWeight<W> &_w1,
                       const VectorWeight<W> &_w2) {
  if (_w1.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER and _w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) return _w1.type_ != _w2.type_;
  const VectorWeight<W> &w1 = __VECTOR_WEIGHT_ASSIGN(_w1, _w2.Length());
  const VectorWeight<W> &w2 = __VECTOR_WEIGHT_ASSIGN(_w2, _w1.Length());

  if (w1.Length() != w2.Length()) {
    LOG(ERROR) << " comparing vectors of different length ";
    return true;
  }

  bool not_equal = false;
  for (size_t i = 0; (i < w1.Length()) && !not_equal; ++i)
    not_equal = not_equal || (w1.Value(i) != w2.Value(i));
  return not_equal;
}

template <class W>
inline bool ApproxEqual(const VectorWeight<W> &_w1,
                        const VectorWeight<W> &_w2,
                        float delta = fst::kDelta) {
  if (_w1.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER and _w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) return _w1.type_ == _w2.type_;
  const VectorWeight<W> &w1 = __VECTOR_WEIGHT_ASSIGN(_w1, _w2.Length());
  const VectorWeight<W> &w2 = __VECTOR_WEIGHT_ASSIGN(_w2, _w1.Length());

  if (w1.Length() != w2.Length()) {
    LOG(ERROR) << " comparing vectors of different length ";
    return false;
  }

  bool approx_equal = true;
  for (size_t i = 0; i < w1.Length(); ++i)
    approx_equal = approx_equal &&
        ApproxEqual(w1.Value(i), w2.Value(i), delta);
  return approx_equal;
}

template <class W>
inline std::ostream &operator<<(std::ostream &strm, const VectorWeight<W> &w) {
  CHECK(FLAGS_fst_weight_separator.size() == 1);
  char separator = FLAGS_fst_weight_separator[0];
  bool write_parens = false;
  if (!FLAGS_fst_weight_parentheses.empty()) {
    CHECK(FLAGS_fst_weight_parentheses.size() == 2);
    write_parens = true;
  }

  if (write_parens)
    strm << FLAGS_fst_weight_parentheses[0];
  for (size_t i = 0; i < w.Length(); ++i) {
    if(i)
      strm << separator;
    strm << w.Value(i);
  }
  if (write_parens)
    strm << FLAGS_fst_weight_parentheses[1];

  return strm;
}

typedef VectorWeight<fst::LogWeight> LogLinearWeight;
typedef fst::ArcTpl<LogLinearWeight> LogLinearArc;
typedef fst::VectorFst<LogLinearArc> LogLinearFst;

}  // namespace openlat

#endif /* openlat_VECTOR_WEIGHT_H_ */
