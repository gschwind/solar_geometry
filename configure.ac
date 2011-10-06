#                     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([solar_geometry], [1.0.3])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/solar_geometry.c])
#AC_CONFIG_HEADER([config.h])

# --enable-m64 compiles in 64 bits. Default is 32 bits as Helioclim does not support 64 bits (yet).
# Order is not arbitrary. Keep this block before AC_PROG_CC/AC_PROG_CXX.
AC_ARG_ENABLE(m64,
[  --enable-m64            compiles in 64 bits],[CFLAGS="$CFLAGS -m64 -mfpmath=387";],[CFLAGS="$CFLAGS -m32";])

# --enable-debug compiles with debugging info. Default is -O2 optimization.
# Order is not arbitrary. Keep this block before AC_PROG_CC/AC_PROG_CXX.
AC_ARG_ENABLE(debug,
[  --enable-debug          compiles with debugging info],[CFLAGS="$CFLAGS -g";],[CFLAGS="$CFLAGS -O2";])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_CONST

# Checks for library functions.
#AC_CHECK_FUNCS([sqrt])

# Turn on gcc checks
if test "x$GCC" = "xyes"; then
  CFLAGS="-fno-strict-aliasing -Wall -Wextra -Wno-unused -Wno-sign-compare -Wstrict-prototypes -Wformat -Wswitch-enum $CFLAGS"
fi

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  solar_geometry.pc
])
AC_OUTPUT
