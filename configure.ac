AC_INIT([BMN],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#silent automake
AM_SILENT_RULES([yes])

#check c++11
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11(noext,mandatory)

#check eigen
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([Eigen/Dense],,AC_MSG_ERROR([Error! Eigen required]))
AC_LANG_POP([C++])

AC_CONFIG_FILES(Makefile)

AC_OUTPUT
