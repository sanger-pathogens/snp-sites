AC_INIT([snp-sites], [1])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_HOST

case $host_os in
	*linux*)
		HOST_OS=linux
		;;
	*freebsd*)
		HOST_OS=freebsd
		;;
	*) ;;
esac

AM_CONDITIONAL([HOST_LINUX],[test x$HOST_OS = xlinux])
AM_CONDITIONAL([HOST_FREEBSD],[test x$HOST_OS = xfreebsd])

AC_CHECK_HEADERS([zlib.h math.h])

AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CXX

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
