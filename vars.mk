# Modified by Vincent Loechner, may 2000.

# Copyright (c) 2000 The Regents of the University of California.
# All rights reserved.
# Permission is hereby granted, without written agreement and without
# license or royalty fees, to use, copy, modify, and distribute this
# software and its documentation for any purpose, provided that the above
# copyright notice and the following two paragraphs appear in all copies
# of this software.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
# ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
# THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
# THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
# PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
# CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
# ENHANCEMENTS, OR MODIFICATIONS.
#
#                                        PT_COPYRIGHT_VERSION_2
#                                        COPYRIGHTENDKEY
#
# Version identification:
# $Id: vars.mk,v 1.3 2001/11/19 13:24:10 loechner Exp $
# Date of creation: 7/31/96
# Author: Bart Kienhuis

VERSION = 5.0

# NOTE: Don't edit this file if it is called vars.mk, instead
# edit vars.mk.in, which is read by configure

# Default top-level directory.
prefix =	/usr

# Usually the same as prefix. 
# exec_prefix is part of the autoconf standard.
exec_prefix =	${prefix}

# Source directory we are building from.
srcdir =	.

# Directory in which to install scripts
BIN_INSTALL_DIR =	$(exec_prefix)/bin

# c compiler flags and defines
CC	             	= gcc
CFLAGS                  = -g -O2
EXTRA_FLAGS             = 
DEFINES         	=  -DSTDC_HEADERS=1 -DHAVE_LIMITS_H=1 -DHAVE_UNISTD_H=1 -DSIZEOF_INT=4 -DSIZEOF_LONG_INT=4 -DSIZEOF_LONG_LONG_INT=8  

# Linker flags and defines, both static and shared lib
LD 	        	= gcc
LDFLAGS                 = 

# Additional libraries that need to be linked
LIBS                    = 

# Additional tools needed when making libraries
LN_S			= ln -s
RANLIB			= ranlib

SHEXT  			= so

# defines needed for arithmetic lib
INT_AFLAGS = -DLINEAR_VALUE_IS_INT
LONG_AFLAGS = -DLINEAR_VALUE_PROTECT_MULTIPLY 		-DLINEAR_VALUE_ASSUME_SOFTWARE_IDIV -DLINEAR_VALUE_IS_LONGLONG
GMP_AFLAGS = -DGNUMP
INT_BITS = 32
LONG_BITS = 64
GMP_BITS = gmp

# Library type to construct
INSTALL_LIB = 

# Commands used to install scripts and data
INSTALL =		/usr/bin/install -c
INSTALL_PROGRAM =	${INSTALL}
INSTALL_DATA =		${INSTALL} -m 644

## GNU-MP stuff
EXTRA_INCLUDES=
EXTRA_LIBS=

# Platform specific variables
OSTYPE	= linux-gnu
HOST    = pc
BUILD   = i686


EXEC_EXTRA_SUFFIX = 

## make install puts everything here: relays on --prefix 
INSTALLDIR = /usr
BINDIR = $(INSTALLDIR)/bin
LIBDIR = $(INSTALLDIR)/lib
INCLUDEDIR = $(INSTALLDIR)/include
MANDIR = $(INSTALLDIR)/man
DOCSDIR = $(INSTALLDIR)/doc/packages/polylib-$(VERSION)

# When compiling the tests, we need to link additional libraries
# include polylib
EXEC_EXTRA_LIBS= -L$(OBJ_DIR) -lpolylib$(BITS) $(LIBS)
SHAREDLIB_FLAG          = -shared
LDCONFIG = ldconfig

LIBS_TO_BUILD = 64
EXEC_TO_BUILD = 64
