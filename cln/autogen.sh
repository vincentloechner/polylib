#!/bin/sh

# MacOS:
if test `uname` = "Darwin"
then
  glibtoolize -c --force
else
  libtoolize -c --force
fi

aclocal -I m4
autoheader
automake -a -c --foreign
autoconf
