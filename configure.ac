AC_INIT([snp-sites], [1])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_HOST

PKG_CHECK_MODULES([CHECK],[check >= 0.8.2],[have_check="yes"],
  AC_MSG_WARN(['Check' unit testing framework not found. It would be impossible to run unit tests!"])
  [have_check="no"])

AC_CHECK_HEADERS([zlib.h math.h])
AC_CHECK_LIB(m, floor)
AX_CHECK_ZLIB

AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CXX

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
