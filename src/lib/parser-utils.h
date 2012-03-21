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
 * parser-utils.h
 *
 *  Created on: 21/03/2012
 *      Author: valabau
 */

#ifndef openlat_PARSER_UTILS_H_
#define openlat_PARSER_UTILS_H_

namespace openlat {


template <typename Context>
inline int istream_input(Context &ctx, char *buf, int max_size) {
  if (not ctx.is.eof()) {
    ctx.is.read(buf, max_size);
    return ctx.is.gcount();
  }
  else return 0;
}

// copy text with optional remove quotes
template <typename Char>
Char *parser_txtcpy(const Char *txt, bool remove_quotes) {
  Char *ptr;
  int len = (remove_quotes)? strlen(txt)-2:strlen(txt);
  ptr = new Char[len+1];
  strncpy(ptr,txt + ((remove_quotes)?1:0), len);
  ptr[len] = '\0';
  return ptr;
}

template <class Loc>
void parser_error(const Loc& l, const std::string& m) {
  std::cerr << "parsing error:" << l << ": " << m << "\n";
}




}


#endif /* openlat_PARSER_UTILS_H_ */
