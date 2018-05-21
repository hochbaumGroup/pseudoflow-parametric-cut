from setuptools import setup, Extension
import unittest


# Fixing needed for building shared library under windows.
from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext


class build_ext_new(build_ext):

    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypes)
        return build_ext.build_extension(self, ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return build_ext.get_export_symbols(self, ext)

    def get_ext_filename(self, ext_name):
        if self._ctypes:
            return ext_name.split('.')[-1] + '.so'
        return build_ext.get_ext_filename(self, ext_name)


class CTypes(Extension):
    pass


def test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests')
    return test_suite


def readme():
    with open('README.md') as f:
        return f.read()


def license():
    with open('LICENSE.md') as f:
        return f.read()


extensions = [
    CTypes(
        "pseudoflow.libhpf",
        ["pseudoflow/core/libhpf.c"],
        depends=["pseudoflow/core/libhpf.h"],
        export_symbols=["hpf_solve", 'libfree'],
        # include_dirs=["pseudoflow/core"],
        language='c99',
        extra_compile_args=['-std=c99', '-O3']
    )
]


setup(name='pseudoflow',
      version='0.1dev',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 2.7',
      ],
      description='Pseudoflow algorithm for computing ' +
      '(parametric) minimum cuts',
      keywords=['minimum cut', 'network flow', 'parametric']
      url='https://github.com/quic0/pseudoflow',
      author='Quico Spaen',
      author_email='qspaen@berkeley.edu',
      license=license(),
      long_description=readme(),
      packages=['pseudoflow', 'pseudoflow/python'],
      install_requires=['networkx', 'six'],
      test_suite='setup.test_suite',
      ext_modules=extensions,
      cmdclass={'build_ext': build_ext_new},
      zip_safe=False)
