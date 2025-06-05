from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import pybind11

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path"""
    def __str__(self):
        return pybind11.get_include()

ext_modules = [
    Extension(
        'UF23Fieldpy',                              # name of the module
        ['UF23Fieldpy/LoopcalcUF23Field.cxx', 'UF23Fieldpy/UF23Field.cc'],  # your source files
        include_dirs=[
            get_pybind_include(),
            # add other include dirs, e.g. current directory, or where UF23Field.h lives
            '.',
        ],
        language='c++',
        extra_compile_args=['-std=c++17'],  # match your C++ standard
    ),
]


setup(
    name='UF23Fieldpy',
    version='1.0.0',
    author='Regnier Mathias',
    description='Pybind11 wrapper for UF23Field',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
    install_requires=[],
)
