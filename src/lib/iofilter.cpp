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
 * file.cpp
 *
 *  Created on: 01/03/2012
 *      Author: valabau
 */

#include <cassert>
#include <fstream>
#include <cstring>
#include <openlat/config.h>
#include <openlat/iofilter.h>

#ifdef HAVE_BOOST_IOSTREAMS
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

using namespace boost::iostreams;
#endif


namespace openlat {


typedef enum { STDIO, GZIP, BZIP2, NO_FILTER } file_filter_t;
static const char *stdio_ext[] = { "-", 0 };
static const char *gzip_ext[] = { ".gz", ".gzip", 0 };
static const char *bzip2_ext[] = { ".bz", ".bz2", ".bzip", ".bzip2", 0 };
static const char **filter_ext[] = { stdio_ext, gzip_ext, bzip2_ext, 0 };


file_filter_t guess_filter(const char * fname) {
  if (strcmp(fname, "-") == 0) return STDIO;
  else {
    const int len = strlen(fname);
    for (int i = 1; filter_ext[i] != 0; i++) {
      for (int j = 0; filter_ext[i][j] != 0; j++) {
        int elen = strlen(filter_ext[i][j]);
        if (len < elen) continue;
        if (strcmp(fname + len - elen, filter_ext[i][j]) == 0) return static_cast<file_filter_t>(i);
      }
    }
  }
  return NO_FILTER;
}


void ifilter::open(const char* name, std::ios_base::openmode mode_) {
  file_filter_t filter_type = guess_filter(name);

  if (filter_type == STDIO) {
    rdbuf(std::cin.rdbuf());
  }
  else {
    bool is_binary = (filter_type == GZIP || filter_type == BZIP2);
    std::ios_base::openmode mode = mode_ | (is_binary?std::ios_base::binary:static_cast<std::ios_base::openmode>(0));
    file_ = new std::ifstream(name, mode);

#ifdef HAVE_BOOST_IOSTREAMS
    filtering_streambuf<input> *filter  = new filtering_streambuf<input>();

    if (filter_type == GZIP) {
      filter->push(gzip_decompressor());
    }
    else if (filter_type == BZIP2) {
      filter->push(bzip2_decompressor());
    }

    filter->push(*file_);

    buf_ = filter;
    std::istream::rdbuf(filter);
#else
    if (filter_type == GZIP || filter_type == BZIP2) {
      std::cerr << "Decompression of gzip and bzip2 files is disabled. To enable it, install the iostreams boost library and recompile\n";
      assert(not (filter_type == GZIP || filter_type == BZIP2));
    }

    std::istream::rdbuf(file_->rdbuf());
#endif
  }
}


void ofilter::open(const char* name, std::ios_base::openmode mode_) {
  file_filter_t filter_type = guess_filter(name);

  if (filter_type == STDIO) {
    rdbuf(std::cout.rdbuf());
  }
  else {
    bool is_binary = (filter_type == GZIP || filter_type == BZIP2);
    std::ios_base::openmode mode = mode_ | (is_binary?std::ios_base::binary:static_cast<std::ios_base::openmode>(0));
    file_ = new std::ofstream(name, mode);

#ifdef HAVE_BOOST_IOSTREAMS
    filtering_streambuf<output> *filter  = new filtering_streambuf<output>();

    if (filter_type == GZIP) {
      filter->push(gzip_compressor());
    }
    else if (filter_type == BZIP2) {
      filter->push(bzip2_compressor());
    }

    filter->push(*file_);

    buf_ = filter;
    std::ostream::rdbuf(filter);
#else
    if (filter_type == GZIP || filter_type == BZIP2) {
      std::cerr << "Decompression of gzip and bzip2 files is disabled. To enable it, install the iostreams boost library and recompile\n";
      assert(not (filter_type == GZIP || filter_type == BZIP2));
    }

    std::ostream::rdbuf(file_->rdbuf());
#endif
  }
}

}

