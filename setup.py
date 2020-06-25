#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
from subprocess import check_output
from collections import defaultdict
import numpy
import os, sys

def pkgconfig(*packages, **kw):
 flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
 ret = defaultdict(list)
 ret.update(kw)
 for token in check_output(["pkg-config", "--libs", "--cflags"] + list(packages)).decode('UTF-8').split():
  if token[:2] in flag_map:
   ret[flag_map.get(token[:2])].append(token[2:])
  else:
   ret['extra_link_args'].append(token)
 for k, v in kw.items(): # remove duplicated
  ret[k] = list(set(ret[k]))
 return ret

params = pkgconfig('solar_geometry')
params['extra_compile_args'] = ['-std=c++11']
params["include_dirs"] += [numpy.get_include(), "third-parties/python-bind-helper", "src"]

# This section try to find the required
# static libraries
if "--static" in sys.argv:
    sys.argv.remove("--static")
    new_libraries = []
    for l in params["libraries"]:
        static_lib_found = False
        for p in params["library_dirs"]:
            f = "%s/lib%s.a"%(p,l)
            if os.path.exists(f):
                params["extra_objects"].append(f)
                static_lib_found = True
                break
        if not static_lib_found:
            new_libraries.append(l)

    params["libraries"] = new_libraries

params["libraries"] = new_libraries


solar_geometry_module = Extension('solar_geometry',
    sources = ['python/pysolar_geometry.cxx'],
    **params)

setup(name = 'solar_geometry',
      version = '1.0',
      description = 'Solar Geometry library',
      ext_modules = [solar_geometry_module])

