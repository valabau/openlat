#!/bin/sh

echo "Generating configure files... may take a while."

mkdir -p m4
autoreconf --install --force -I m4 -v && \
  echo "Preparing was successful if there was no error messages above." && \
  echo "Now type:" && \
  echo "  ./configure && make"  && \
  echo "Run './configure --help' for more information"
