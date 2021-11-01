#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_thermo',
                           sources=['thermoInterface_wrap.C', 'thermoInterface.C'],
                           include_dirs=['/home/zhy/new_drive/OpenFOAM-6/user/platforms/linux64GccDPInt32Opt/python'],
                           library_dirs=['/home/zhy/new_drive/OpenFOAM-6/user/platforms/linux64GccDPInt32Opt/python'],
                           libraries = ["_thermo"]
                           )

setup (name = 'thermo',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [example_module],
       py_modules = ["example"],
       )