AC_INIT([polylib],[5.22.8])
AM_INIT_AUTOMAKE([gnu subdir-objects])
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX
LT_INIT

AC_LANG(C++)

AC_PATH_CLN(1.1.0, [
  LIBS="$LIBS $CLN_LIBS"
  CPPFLAGS="$CPPFLAGS $CLN_CPPFLAGS"
], AC_MSG_ERROR([No suitable installed version of CLN could be found.]))

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES([Makefile polylib/Makefile])
AC_OUTPUT
