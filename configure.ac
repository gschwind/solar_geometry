# -*- coding: utf-8 -*-
#
#                     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([solar_geometry], [0.0.1])

AC_CONFIG_MACRO_DIR([m4])

LT_PREREQ([2.4.2])
LT_INIT()

AM_INIT_AUTOMAKE()
AM_MAINTAINER_MODE()
AC_CONFIG_SRCDIR([src/solar_geometry.cxx])

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL

AX_CXX_COMPILE_STDCXX_11(noext, mandatory)

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  solar_geometry.pc
])

AC_OUTPUT
