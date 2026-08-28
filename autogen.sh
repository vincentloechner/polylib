#!/bin/sh

# MacOS:
if test `uname` = "Darwin"
then
  glibtoolize -c --force
fi

# libtoolize -c --force
# aclocal -I m4
# autoheader
# automake -a -c --foreign
# autoconf

autoreconf -i
(cd cln; ./autogen.sh)
