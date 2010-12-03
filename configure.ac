#                     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([solar_geometry], [1.0.3])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/solar_geometry.c])
#AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_CONST

# Checks for library functions.
#AC_CHECK_FUNCS([sqrt])

# Default gcc options
if test "x$GCC" = "xyes"; then
  CFLAGS="-fno-strict-aliasing -Wall -Wextra -Wno-unused -Wno-sign-compare -Wstrict-prototypes -Wformat -Wswitch-enum -Werror $CFLAGS"
fi

# --enable-debug compiles with debugging info
AC_ARG_ENABLE(debug,
[  --enable-debug          compiles with debugging info],
[if test x$enableval = xyes; then
	CFLAGS="$CFLAGS -g -O0"
fi])

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  solar_geometry.pc
])
AC_OUTPUT
