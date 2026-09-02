#!/bin/sh

# MacOS:
if test `uname` = "Darwin"
then
  glibtoolize -c --force
fi

# libtoolize -c --force
# aclocal -I m4
# autoheader
# automake -a -c --foreign --add-missing
# autoconf

autoreconf -i

echo "Autogen in cln/:"
(cd cln; ./autogen.sh)
