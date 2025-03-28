# Configure paths for the CLN library
# Richard Kreckel 12/4/2000
# borrowed from Christian Bauer
# stolen from Sam Lantinga
# stolen from Manish Singh
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl AC_PATH_CLN([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for installed CLN library, and define CLN_CPPFLAGS and CLN_LIBS
dnl
AC_DEFUN([AC_PATH_CLN],
[dnl 
dnl Get the cppflags and libraries from the cln-config script
dnl
AC_ARG_WITH(cln-prefix,[  --with-cln-prefix=PFX   Prefix where CLN is installed (optional)],
            cln_config_prefix="$withval", cln_config_prefix="")
AC_ARG_WITH(cln-exec-prefix,[  --with-cln-exec-prefix=PFX Exec prefix where CLN is installed (optional)],
            cln_config_exec_prefix="$withval", cln_config_exec_prefix="")
AC_ARG_ENABLE(clntest, [  --disable-clntest       Do not try to compile and run a test CLN program],
              , enable_clntest=yes)

if test x$cln_config_exec_prefix != x ; then
    cln_config_args="$cln_config_args --exec-prefix=$cln_config_exec_prefix"
    if test x${CLN_CONFIG+set} != xset ; then
        CLN_CONFIG=$cln_config_exec_prefix/bin/cln-config
    fi
fi
if test x$cln_config_prefix != x ; then
    cln_config_args="$cln_config_args --prefix=$cln_config_prefix"
    if test x${CLN_CONFIG+set} != xset ; then
        CLN_CONFIG=$cln_config_prefix/bin/cln-config
    fi
fi

AC_PATH_PROG(CLN_CONFIG, cln-config, no)
cln_min_version=ifelse([$1], ,1.1.0,$1)
AC_MSG_CHECKING(for CLN - version >= $cln_min_version)
if test "$CLN_CONFIG" = "no" ; then
    AC_MSG_RESULT(no)
    echo "*** The cln-config script installed by CLN could not be found"
    echo "*** If CLN was installed in PREFIX, make sure PREFIX/bin is in"
    echo "*** your path, or set the CLN_CONFIG environment variable to the"
    echo "*** full path to cln-config."
    ifelse([$3], , :, [$3])
else
dnl Parse required version and the result of cln-config.
    cln_min_major_version=`echo $cln_min_version | \
            sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    cln_min_minor_version=`echo $cln_min_version | \
            sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    cln_min_micro_version=`echo $cln_min_version | \
            sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    CLN_CPPFLAGS=`$CLN_CONFIG $cln_config_args --cppflags`
    CLN_LIBS=`$CLN_CONFIG $cln_config_args --libs`
    cln_config_version=`$CLN_CONFIG $cln_config_args --version`
    cln_config_major_version=`echo $cln_config_version | \
            sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    cln_config_minor_version=`echo $cln_config_version | \
            sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    cln_config_micro_version=`echo $cln_config_version | \
            sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
dnl Check if the installed CLN is sufficiently new according to cln-config.
    if test \( $cln_config_major_version -lt $cln_min_major_version \) -o \
            \( $cln_config_major_version -eq $cln_min_major_version -a $cln_config_minor_version -lt $cln_min_minor_version \) -o \
            \( $cln_config_major_version -eq $cln_min_major_version -a $cln_config_minor_version -eq $cln_min_minor_version -a $cln_config_micro_version -lt $cln_min_micro_version \); then
        echo -e "\n*** 'cln-config --version' returned $cln_config_major_version.$cln_config_minor_version.$cln_config_micro_version, but the minimum version"
        echo "*** of CLN required is $cln_min_major_version.$cln_min_minor_version.$cln_min_micro_version. If cln-config is correct, then it is"
        echo "*** best to upgrade to the required version."
        echo "*** If cln-config was wrong, set the environment variable CLN_CONFIG"
        echo "*** to point to the correct copy of cln-config, and remove the file"
        echo "*** config.cache before re-running configure."
        ifelse([$3], , :, [$3])
    else
dnl The versions match so far.  Now do a sanity check: Does the result of cln-config
dnl match the version of the headers and the version built into the library, too?
        no_cln=""
        if test "x$enable_clntest" = "xyes" ; then
            ac_save_CPPFLAGS="$CPPFLAGS"
            ac_save_LIBS="$LIBS"
            CPPFLAGS="$CPPFLAGS $CLN_CPPFLAGS"
            LIBS="$LIBS $CLN_LIBS"
            rm -f conf.clntest
            AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdio.h>
#include <string.h>
#include <cln/version.h>

/* we do not #include <stdlib.h> because autoconf in C++ mode inserts a
   prototype for exit() that conflicts with the one in stdlib.h */
extern "C" int system(const char *);

int main(void)
{
    int major, minor, micro;
    char *tmp_version;

    system("touch conf.clntest");

    if ((CL_VERSION_MAJOR != $cln_config_major_version) ||
        (CL_VERSION_MINOR != $cln_config_minor_version) ||
        (CL_VERSION_PATCHLEVEL != $cln_config_micro_version)) {
        printf("\n*** 'cln-config --version' returned %d.%d.%d, but the header file I found\n", $cln_config_major_version, $cln_config_minor_version, $cln_config_micro_version);
        printf("*** corresponds to %d.%d.%d. This mismatch suggests your installation of CLN\n", CL_VERSION_MAJOR, CL_VERSION_MINOR, CL_VERSION_PATCHLEVEL);
        printf("*** is corrupted or you have specified some wrong -I compiler flags.\n");
        printf("*** Please inquire and consider reinstalling CLN.\n");
        return 1;
    }
    if ((cln::version_major != $cln_config_major_version) ||
        (cln::version_minor != $cln_config_minor_version) ||
        (cln::version_patchlevel != $cln_config_micro_version)) {
        printf("\n*** 'cln-config --version' returned %d.%d.%d, but the library I found\n", $cln_config_major_version, $cln_config_minor_version, $cln_config_micro_version);
        printf("*** corresponds to %d.%d.%d. This mismatch suggests your installation of CLN\n", cln::version_major, cln::version_minor, cln::version_patchlevel);
        printf("*** is corrupted or you have specified some wrong -L compiler flags.\n");
        printf("*** Please inquire and consider reinstalling CLN.\n");
        return 1;
    }
    return 0;
}
]])],[],[no_cln=yes],[echo $ac_n "cross compiling; assumed OK... $ac_c"])
            CPPFLAGS="$ac_save_CPPFLAGS"
            LIBS="$ac_save_LIBS"
        fi
        if test "x$no_cln" = x ; then
            AC_MSG_RESULT([yes, $cln_config_version])
            ifelse([$2], , :, [$2])
        else
            AC_MSG_RESULT(no)
            if test ! -f conf.clntest ; then
                echo "*** Could not run CLN test program, checking why..."
                CPPFLAGS="$CFLAGS $CLN_CPPFLAGS"
                LIBS="$LIBS $CLN_LIBS"
                AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <stdio.h>
#include <cln/version.h>
]], [[ return 0; ]])],[ echo "*** The test program compiled, but did not run. This usually means"
                  echo "*** that the run-time linker is not finding CLN or finding the wrong"
                  echo "*** version of CLN. If it is not finding CLN, you'll need to set your"
                  echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
                  echo "*** to the installed location. Also, make sure you have run ldconfig if that"
                  echo "*** is required on your system."],[ echo "*** The test program failed to compile or link. See the file config.log for the"
                  echo "*** exact error that occured. This usually means CLN was incorrectly installed"
                  echo "*** or that you have moved CLN since it was installed. In the latter case, you"
                  echo "*** may want to edit the cln-config script: $CLN_CONFIG." ])
                CPPFLAGS="$ac_save_CPPFLAGS"
                LIBS="$ac_save_LIBS"
            fi
            CLN_CPPFLAGS=""
            CLN_LIBS=""
            ifelse([$3], , :, [$3])
        fi
    fi
fi
AC_SUBST(CLN_CPPFLAGS)
AC_SUBST(CLN_LIBS)
rm -f conf.clntest
])
