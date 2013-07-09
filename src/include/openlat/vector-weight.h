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

void dump_backtrace();

namespace openlat {

// vector weight
template <class W>
class VectorWeight {
  typedef enum { VECTOR_WEIGHT_ZERO, VECTOR_WEIGHT_ONE, VECTOR_WEIGHT_OTHER, VECTOR_WEIGHT_NONE } vector_weight_t;

  template<class U> friend bool operator== (const VectorWeight<U> &_w1, const VectorWeight<U> &_w2);
  template<class U> friend bool operator!= (const VectorWeight<U> &_w1, const VectorWeight<U> &_w2);
  template<class U> friend bool ApproxEqual(const VectorWeight<U> &_w1, const VectorWeight<U> &_w2, float delta);
  template<class U> friend inline VectorWeight<U> Plus(const VectorWeight<U> &w1, const VectorWeight<U> &w2);
  template<class U> friend inline VectorWeight<U> Times(const VectorWeight<U> &w1, const VectorWeight<U> &w2);
  template<class U> friend inline VectorWeight<U> Divide(const VectorWeight<U> &w1, const VectorWeight<U> &w2, fst::DivideType typ);

public:
  typedef W Weight;
  typedef VectorWeight<typename W::ReverseWeight> ReverseWeight;

  void check() {
    if ((values_.size() == 0 and type_ == VECTOR_WEIGHT_OTHER) or (values_.size() > 0 and type_ != VECTOR_WEIGHT_OTHER)) {
      LOG(ERROR) << "Invalid VectorWeight "<< Length() << "(" << type_ << ") - " << combination_;
      for (size_t i = 0; i < Length(); i++) LOG(ERROR) << "  w[" << i << "] = " << values_[i];
      dump_backtrace();
    }
  }
  VectorWeight(vector_weight_t type = VECTOR_WEIGHT_ONE): values_(0), type_(type) {
    if (type == VECTOR_WEIGHT_OTHER) type = VECTOR_WEIGHT_ONE;
    combination_ = (type == VECTOR_WEIGHT_NONE)?Weight::NoWeight():((type == VECTOR_WEIGHT_ZERO)?Weight::Zero():Weight::One());
    check();
  }

  VectorWeight(const VectorWeight &w): combination_(w.combination_), values_(w.values_), type_(w.type_) {check();}

//  VectorWeight(const W& c, const vector<W> &w): combination_(c), values_(w), type_(VECTOR_WEIGHT_OTHER) {check();}

  template <typename T>
  VectorWeight(const T& c,const vector<T> &w): type_(VECTOR_WEIGHT_OTHER) {
    combination_ = W(c);
    Resize(w.size());
    for (size_t i = 0; i < w.size(); i++) values_[i] = W(w[i]);
    check();
  }

//  template <class Iterator>
//  VectorWeight(const W& c,Iterator begin, Iterator end): combination_(c), values_(begin, end), type_(VECTOR_WEIGHT_OTHER) {check();}

//  VectorWeight(const W& c, const W &w, size_t n): combination_(c), values_(vector<W>(n, w)), type_(VECTOR_WEIGHT_OTHER) {check();}

  static uint64 Properties() {
    return Weight::Properties();
  }

  void Resize(size_t n) {
    type_ = VECTOR_WEIGHT_OTHER;
    if (type_ == VECTOR_WEIGHT_ZERO) {
      values_.clear();
      values_.resize(n, W::Zero());
    }
    else if (type_ == VECTOR_WEIGHT_ONE) {
      values_.clear();
      values_.resize(n, W::One());
    }
    else if (type_ == VECTOR_WEIGHT_NONE) {
      values_.clear();
      values_.resize(n, W::NoWeight());
    }
    else {
      values_.resize(n, W::One());
    }
  }

  size_t Length() const {
    return values_.size();
  }

  void Update(const vector<float> &w) {
    if (type_ != VECTOR_WEIGHT_OTHER) return;

    if (Length() > w.size()) LOG(ERROR) << " updating with weights of different lengths ";

    Resize(w.size());

    combination_ = W::One();
    for (size_t i = 0; i < w.size(); i++) {
      combination_ = Times(combination_, W(values_[i].Value() * w[i]));
    }
  }

  static const VectorWeight<W> &Zero() {
    static const VectorWeight<W> zero(VECTOR_WEIGHT_ZERO);
    return zero;
  }

  static const VectorWeight<W> &One() {
    static const VectorWeight<W> one(VECTOR_WEIGHT_ONE);
    return one;
  }

  static const VectorWeight<W> NoWeight() {
    static const VectorWeight<W> none(VECTOR_WEIGHT_NONE);
    return none;
  }

  static const ReverseWeight &ReverseZero() {
    static const ReverseWeight zero(VECTOR_WEIGHT_ZERO);
    return zero;
  }

  static const ReverseWeight &ReverseOne() {
    static const ReverseWeight one(VECTOR_WEIGHT_ONE);
    return one;
  }

  static const ReverseWeight ReverseNoWeight() {
    static const ReverseWeight none(VECTOR_WEIGHT_NONE);
    return none;
  }


  std::ostream &Write(std::ostream &strm) const {
    combination_.Write(strm);
    for (size_t i = 0; i < Length(); ++i)
      values_[i].Write(strm);
    return strm;
  }

  VectorWeight<W> &operator=(const VectorWeight<W> &w) {
    type_ = w.type_;
    combination_ = w.combination_;
    values_ = w.values_;
    check();
    return *this;
  }

  bool Member() const {
    // Zero and One are already memebers
    if (type_ == VECTOR_WEIGHT_ZERO) return W::Zero().Member();
    else if (type_ == VECTOR_WEIGHT_ONE) return W::One().Member();
    else if (type_ == VECTOR_WEIGHT_NONE) return W::NoWeight().Member();
    else {
      bool member = combination_.Member();
      for (size_t i = 0; i < Length(); ++i)
        member = member && values_[i].Member();
      return member;
    }
  }

  size_t Hash() const {
    if (type_ == VECTOR_WEIGHT_ZERO) return size_t(W::Zero().Hash());
    else if (type_ == VECTOR_WEIGHT_ONE) return size_t(W::One().Hash());
    else if (type_ == VECTOR_WEIGHT_NONE) return size_t(W::NoWeight().Hash());
    else {
      uint64 hash = combination_.Hash();
      for (size_t i = 0; i < Length(); ++i)
        hash = 5 * hash + values_[i].Hash();
      return size_t(hash);
    }
  }

  VectorWeight<W> Quantize(float delta = fst::kDelta) const {
    if (type_ == VECTOR_WEIGHT_OTHER) {
      VectorWeight<W> w;
      w.type_ = type_;
      w.combination_ = combination_.Quantize(delta);
      w.Resize(Length());
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
    if (type_ == VECTOR_WEIGHT_NONE) {
      return ReverseNoWeight();
    }
    else {
      ReverseWeight w;
      w.type_ = type_;
      w.combination_ = combination_.Reverse();
      w.Resize(Length());
      for (size_t i = 0; i < Length(); ++i)
        w.values_[i] = values_[i].Reverse();
      return w;
    }
  }

  const float& Value() const {
    return combination_.Value();
  }

  const float& Value(size_t i) const {
    static const W zero(W::Zero());
    static const W one(W::One());
    static const W none(W::NoWeight());
    if (type_ == VECTOR_WEIGHT_ZERO)      return zero.Value();
    else if (type_ == VECTOR_WEIGHT_ONE)  return one.Value();
    else if (type_ == VECTOR_WEIGHT_NONE) return none.Value();
    else if (i < Length()) {
      return values_[i].Value();
    }
    else {
      LOG(ERROR) << " accessing to vector weight outside boundaries at index "
                 << i << " out of " << Length()
                 << " for type " << type_;
      return one.Value();
    }
  }

  static const string &Type(void) {
    static const string type = "vector_" + W::Type();
    return type;
  }


 protected:

 private:
  W combination_;
  vector<W> values_;
  vector_weight_t type_;

};

template <class W>
inline bool operator==(const VectorWeight<W> &w1,
                       const VectorWeight<W> &w2) {
  if (w1.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER and w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) return w1.type_ == w2.type_;
  if (w1.Value() != w2.Value()) return false;
  else {
    if (w1.type_ == w2.type_ and w1.Length() != w2.Length()) {
      LOG(ERROR) << " in operator== comparing vectors of different lengths " << w1.Length() << "(" << w1.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
      return false;
    }

    for (size_t i = 0; i < w1.Length(); ++i)
      if (w1.Value(i) != w2.Value(i)) return false;
    return true;
  }
}

template <class W>
inline bool operator!=(const VectorWeight<W> &w1,
                       const VectorWeight<W> &w2) {
  if (w1.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER and w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) return w1.type_ != w2.type_;
  if (w1.Value() != w2.Value()) return true;
  else {
    if (w1.type_ == w2.type_ and w1.Length() != w2.Length()) {
      LOG(ERROR) << " in operator!= comparing vectors of different lengths " << w1.Length() << "(" << w1.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
      return true;
    }

    for (size_t i = 0; i < w1.Length(); ++i)
      if (w1.Value(i) != w2.Value(i)) return true;
    return false;
  }
}

template <class W>
inline bool ApproxEqual(const VectorWeight<W> &w1,
                        const VectorWeight<W> &w2,
                        float delta = fst::kDelta) {

  if (w1.type_ == w2.type_ and w1.Length() != w2.Length()) {
    LOG(ERROR) << " in ApproxEqual comparing vectors of different lengths " << w1.Length() << "(" << w1.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
    return false;
  }

  if (not ApproxEqual(w1.combination_, w2.combination_, delta)) return false;

  for (size_t i = 0; i < w1.Length(); ++i)
    if (not ApproxEqual(w1.values_[i], w2.values_[i], delta)) return false;
  return true;
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
  strm << w.Value();
  for (size_t i = 0; i < w.Length(); ++i) {
    strm << separator;
    strm << w.Value(i);
  }
  if (write_parens)
    strm << FLAGS_fst_weight_parentheses[1];

  return strm;
}

template <class W>
inline VectorWeight<W> Plus(const VectorWeight<W> &w1,
                             const VectorWeight<W> &w2) {
  VectorWeight<W> res;
  res.combination_ = Plus(w1.combination_, w2.combination_);
  if (w1.combination_.Value() < w2.combination_.Value()) {
    res.values_ = w1.values_;
    res.type_ = w1.type_;
  }
  else {
    res.values_ = w2.values_;
    res.type_ = w2.type_;
  }
  return res;
}

template <class W>
inline VectorWeight<W> Times(const VectorWeight<W> &w1,
                             const VectorWeight<W> &w2) {
  switch (w1.type_) {
    case VectorWeight<W>::VECTOR_WEIGHT_ONE: return w2;
    case VectorWeight<W>::VECTOR_WEIGHT_ZERO:;
    case VectorWeight<W>::VECTOR_WEIGHT_NONE: return w1;
    default: {
      switch (w2.type_) {
        case VectorWeight<W>::VECTOR_WEIGHT_ONE: return w1;
        case VectorWeight<W>::VECTOR_WEIGHT_ZERO:;
        case VectorWeight<W>::VECTOR_WEIGHT_NONE: return w2;
      }
    }
  }

  // else both cases are OTHER
  if (w1.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER or w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) {
    LOG(ERROR) << " in Times comparing vectors of wrong types " << w1.Length() << "(" << w1.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
    return VectorWeight<W>::NoWeight();
  }
  else if (w1.Length() != w2.Length()) {
    LOG(ERROR) << " in Times comparing vectors of different lengths " << w1.Length() << "(" << w1.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
    return VectorWeight<W>::NoWeight();
  }
  else {
    VectorWeight<W> res;
    res.type_ = VectorWeight<W>::VECTOR_WEIGHT_OTHER;
    res.combination_ = Times(w1.combination_, w2.combination_);
    res.values_.resize(w1.Length());
    for (size_t i = 0; i < res.Length(); i++) {
      res.values_[i] = Times(w1.values_[i], w2.values_[i]);
    }
    return res;
  }
}

template <class W>
inline VectorWeight<W> Divide(const VectorWeight<W> &w1,
                              const VectorWeight<W> &w2,
                                    fst::DivideType typ = fst::DIVIDE_ANY) {

  switch (w2.type_) {
    case VectorWeight<W>::VECTOR_WEIGHT_ONE: return w1;
    case VectorWeight<W>::VECTOR_WEIGHT_ZERO:;
    case VectorWeight<W>::VECTOR_WEIGHT_NONE: return VectorWeight<W>::NoWeight();
    default: {
      switch (w1.type_) {
        case VectorWeight<W>::VECTOR_WEIGHT_ZERO: return VectorWeight<W>::Zero();
        case VectorWeight<W>::VECTOR_WEIGHT_NONE: return w1;
        //case VectorWeight<W>::VECTOR_WEIGHT_ONE: w1.Resize(w2.Length());
      }
    }
  }

  VectorWeight<W> res(w1);
  res.Resize(w2.Length());

  // else both cases are OTHER
  if (res.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER or w2.type_ != VectorWeight<W>::VECTOR_WEIGHT_OTHER) {
    LOG(ERROR) << " in Times comparing vectors of wrong types " << res.Length() << "(" << res.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
    return VectorWeight<W>::NoWeight();
  }
  else if (res.Length() != w2.Length()) {
    LOG(ERROR) << " in Times comparing vectors of different lengths " << res.Length() << "(" << res.type_ << ") != " << w2.Length() << "(" << w2.type_ << ")";
    return VectorWeight<W>::NoWeight();
  }
  else {
    res.type_ = VectorWeight<W>::VECTOR_WEIGHT_OTHER;
    res.combination_ = Divide(res.combination_, w2.combination_);
    res.values_.resize(res.Length());
    for (size_t i = 0; i < res.Length(); i++) {
      res.values_[i] = Divide(res.values_[i], w2.values_[i]);
    }
    return res;
  }
}

typedef VectorWeight<fst::TropicalWeight> LogLinearWeight;
typedef fst::ArcTpl<LogLinearWeight> LogLinearArc;
typedef fst::VectorFst<LogLinearArc> LogLinearFst;

}  // namespace openlat

#endif /* openlat_VECTOR_WEIGHT_H_ */
