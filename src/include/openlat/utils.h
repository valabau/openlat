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
 * utils.h
 *
 *  Created on: 17/02/2012
 *      Author: valabau
 */

#ifndef openlat_UTILS_H_
#define openlat_UTILS_H_

#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <cerrno>
#include <cmath>

namespace openlat {


#define print_stack_trace
#define _WRITE(msg) std::cerr << msg

#define _WRITE_MSG(file, line, func_name, expstr, msg) \
do { \
  _WRITE(file << " (" << line << "): "); \
  if (std::string(func_name).length() > 0) { \
    _WRITE(" in function '" << func_name << "': "); \
  } \
  if (std::string(expstr).length() > 0) { \
    _WRITE(" expression (" << expstr << ") failed: "); \
  } \
  if (std::string(#msg).length() > 0) _WRITE(msg); \
  _WRITE("\n"); \
} while(0)


#ifdef  NDEBUG
# define assert_bt(expr, ...)           (__ASSERT_VOID_CAST (0))
#else
# define assert_bt(expr, ...) do { \
    if (!(expr)) { \
      print_stack_trace(2); \
      _WRITE_MSG(__FILE__ , __LINE__, __FUNCTION__, #expr, ##__VA_ARGS__); \
      if (errno > 0) { \
        _WRITE("Critical error: " << strerror(errno) << "\n"); \
      } \
      assert(expr); \
    } \
  } while(0)
#endif


#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif


template<class T>
inline T add_log(T a, T b) {
 /* Make a the max value and b the min value */
 if (b > a) { T tmp = a; a = b; b = tmp; }
 /* If min is zero, return max */
 if (b == -std::numeric_limits<T>::infinity()) { return a; }
 /* Robust sum */
 return a + log(1.0 + exp(b-a));
}


template<class T>
inline T sub_log(T a, T b) {
  /* Make a the max value and b the min value */
  if (b > a) { T tmp = a; a = b; b = tmp; }
  /* If min is zero, return max */
  if (b == -std::numeric_limits<T>::infinity()) { return a; }
  /* Robust sum */
  return a + log(1.0 - exp(b-a));
}

template<typename T>
T change_base(T value, T base) {
  return value * log(base);
}

template<class T> unsigned int edit_distance(const T& s1, const T& s2) {
  const size_t len1 = s1.size(), len2 = s2.size();
  std::vector<vector<unsigned int> > d(len1 + 1, std::vector<unsigned int> (len2 + 1));

  for (size_t i = 1; i <= len1; ++i) d[i][0] = i;
  for (size_t i = 1; i <= len2; ++i) d[0][i] = i;

  for (size_t i = 1; i <= len1; ++i)
    for (size_t j = 1; j <= len2; ++j)
      d[i][j] = std::min(std::min(d[i - 1][j] + 1, d[i][j - 1] + 1), d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1));
  return d[len1][len2];
}

template<class T> unsigned int hamming_distance(const T& s1, const T& s2) {
  assert_bt(s1.size() == s2.size(), "Sizes do not match");

  size_t distance = 0;
  for (size_t i = 0; i < s1.size(); i++) {
    if (s1[i] != s2[i]) distance++;
  }
  return distance;
}

template <typename T>
T convert_string(const string& text) {
  T ret;
  std::istringstream ss(text, std::istringstream::in);
  ss >> ret;
  if (ss.fail()) {
    if (text == "inf") return std::numeric_limits<T>::infinity();
    else if (text == "-inf") return -std::numeric_limits<T>::infinity();
    else if (text == "nan") return -std::numeric_limits<T>::quiet_NaN();
    assert_bt(false, "Unknown error");
  }
  return ret;
}

//template <>
//size_t convert_string<size_t>(const string& text) {
//  size_t ret;
//  std::istringstream ss(text, std::istringstream::in);
//  ss >> ret;
//  if (ss.fail()) assert_bt(false, "Unknown error");
//  return ret;
//}

template <typename T>
string to_string(const T& value) {
  std::ostringstream ss(std::ostringstream::out);
  ss << value;
  return ss.str();
}

template <typename T>
void tokenize(const T& str,
              vector<T>& tokens,
              const T& delimiters = T(" "))
{
  if (delimiters == T("")) {
    tokens.reserve(str.size());
    for (typename T::const_iterator it = str.begin(); it != str.end(); ++it) {
      T tok;
      tok.push_back(*it);
      tokens.push_back(tok);
    }
  }
  else {
    tokens.clear();

    // Skip delimiters at beginning.
    typename T::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    typename T::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (T::npos != pos || T::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
  }
}



/** compute the laps of time in ms between two clocks
 * @param clock1 first clock
 * @param clock2 second clock
 * @return return the lapse of time in ms between clock1 and clock2
 */
template <typename T>
double diffclock(T clock1, T clock2) {
  double diffticks = static_cast<double>(clock1) - static_cast<double>(clock2);
  double diffms = (diffticks * 1000) / CLOCKS_PER_SEC;
  return diffms;
}

template <typename W>
float to_float(const W& weight) {
  return exp(-weight.Value());
}

}

#endif /* openlat_UTILS_H_ */
