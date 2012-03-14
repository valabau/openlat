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
 * file.h
 *
 *  Created on: 24/02/2012
 *      Author: valabau
 */

#ifndef openlat_IOFILTER_H_
#define openlat_IOFILTER_H_

#include <iostream>
#include <fstream>

namespace openlat {

  class ifilter : public std::istream {
  public:
    ifilter(): std::istream(0), file_(0), buf_(0) {};
    ifilter(const char* name, std::ios_base::openmode mode = ios_base::in): std::istream(0), file_(0), buf_(0) { open(name, mode); };
    virtual ~ifilter() { close(); }
    void open(const char* name, std::ios_base::openmode mode = ios_base::in);
    void close() { if (buf_ != 0) { delete buf_; buf_ = 0; }; if (file_ != 0) { file_->close(); delete file_; file_ = 0; }  }
  private:
    std::ifstream *file_;
    std::streambuf *buf_;
  };

  class ofilter : public std::ostream {
  public:
    ofilter(): std::ostream(0), file_(0), buf_(0) {}
    ofilter(const char* name, std::ios_base::openmode mode = ios_base::out): std::ostream(0), file_(0), buf_(0) { open(name, mode); }
    virtual ~ofilter() { close(); }
    void open(const char* name, std::ios_base::openmode mode = ios_base::out);
    void close() { if (buf_ != 0) { delete buf_; buf_ = 0; }; if (file_ != 0) { file_->close(); delete file_; file_ = 0; }; }
    ofilter& flush () { std::ostream::flush(); if (file_ != 0) file_->flush(); else std::cout.flush(); return *this;}
  private:
    std::ofstream *file_;
    std::streambuf *buf_;
  };

  template <typename Context>
  inline int istream_input(Context &ctx, char *buf, int max_size) {
    if (not ctx.is.eof()) {
      ctx.is.read(buf, max_size);
      return ctx.is.gcount();
    }
    else return 0;
  }
}

#endif /* openlat_IOFILTER_H_ */
