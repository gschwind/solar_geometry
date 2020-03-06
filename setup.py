#!/usr/bin/env python3

#from distutils.core import setup, Extension
import subprocess

from distutils.core import setup, Extension
import numpy as np

def pkgconfig(*packages, **kw):
 flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
 for token in subprocess.check_output(["pkg-config", "--libs", "--cflags"] + list(packages)).decode('UTF-8').split():
  if token[:2] in flag_map:
   kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
  else: # throw others to extra_link_args
   kw.setdefault('extra_link_args', []).append(token)
 try:
  # python 2.x
  for k, v in kw.iteritems(): # remove duplicated
   kw[k] = list(set(v))
 except AttributeError:
  # python 3.x
  for k, v in kw.items(): # remove duplicated
   kw[k] = list(set(kw[k]))
 return kw

params = {'extra_link_args': []}
params['extra_compile_args'] = ['-std=c++11', '-ggdb']

if "include_dirs" in params:
 params["include_dirs"] += [np.get_include()]
else:
 params["include_dirs"] = [np.get_include()]

module1 = Extension('solar_geometry', sources = ['python/pysolar_geometry.cxx', 'src/solar_geometry.cxx'], **params)

setup (name = 'solar_geometry',
		version = '1.0',
		description = 'Solar Geometry library',
		ext_modules = [module1])

