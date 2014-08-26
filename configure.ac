# -*- coding: utf-8 -*-
#
#                     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([solar_geometry], [2.0])

LT_PREREQ([2.4.2])
LT_INIT()

AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/solar_geometry.c])
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_HEADER([config.h])

# --enable-debug compiles with debugging info. Default is "-g -O2".
# Order is not arbitrary. Keep this block before AC_PROG_CC/AC_PROG_CXX.
AC_ARG_ENABLE(debug,
 [  --enable-debug          compiles with debugging info],
 [CFLAGS="$CFLAGS -g -O0";],
 [CFLAGS="$CFLAGS -O2 -DNDEBUG";])

# Checks for programs.
AC_PROG_CC
#AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_CONST

# Checks for library functions.
#AC_CHECK_FUNCS([sqrt])

# Use SSE2 floating point instructions in 32 bits mode to get the same results as 64 bits mode.
if test "x$GCC" = "xyes"; then
  CFLAGS="$CFLAGS -mfpmath=sse -msse2"
fi

# Turn on gcc checks
if test "x$GCC" = "xyes"; then
  CFLAGS="-fno-strict-aliasing -Wall -Wextra -Wno-unused -Wno-sign-compare -Wno-comment -Wstrict-prototypes -Wformat -Wswitch-enum $CFLAGS"
fi

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  solar_geometry.pc
])

AC_OUTPUT
