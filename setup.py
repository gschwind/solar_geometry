#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
from subprocess import check_output
import numpy

solar_geometry_module = Extension('solar_geometry',
    sources = ['python/pysolar_geometry.cxx', 'src/solar_geometry.cxx'],
    extra_compile_args = ['-std=c++11', '-ggdb'],
    include_dirs = [numpy.get_include(), "third-parties/python-bind-helper", "src"])

setup(name = 'solar_geometry',
      version = '1.0',
      description = 'Solar Geometry library',
      ext_modules = [solar_geometry_module])

