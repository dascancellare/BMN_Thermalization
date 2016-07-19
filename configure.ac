AC_INIT([BMN],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#silent automake
AM_SILENT_RULES([yes])

AC_PROG_CC

#check c++11
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11(noext,mandatory)

#number of colors
AC_ARG_ENABLE(ncol,
        AS_HELP_STRING([--enable-ncol=n],[Set number of colors [Default: 3]]),
        ncol="${enableval}",
        ncol=3)
if ! printf "%f" "$col" >/dev/null 2>&1; then AC_MSG_ERROR([String $ncol is not a number]); fi
AC_MSG_RESULT([enabling number of colors... $ncol])
AC_DEFINE_UNQUOTED([NCOL],$ncol,[Number of colors])

#check eigen
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([Eigen/Dense],,AC_MSG_ERROR([Error! Eigen required]))
AC_LANG_POP([C++])

AC_CONFIG_FILES(Makefile)

AC_OUTPUT
