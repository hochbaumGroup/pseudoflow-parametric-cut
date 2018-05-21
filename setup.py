from setuptools import setup, Extension, find_packages


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


def readme():
    with open('README.md') as f:
        return f.read()


def license():
    with open('LICENSE.md') as f:
        return f.read()


extensions = [
    CTypes(
        "pseudoflow.libhpf",
        ["src/pseudoflow/core/libhpf.c"],
        depends=["src/pseudoflow/core/libhpf.h"],
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
      keywords=['minimum cut', 'network flow', 'parametric'],
      url='https://github.com/quic0/pseudoflow',
      author='Quico Spaen',
      author_email='qspaen@berkeley.edu',
      license=license(),
      long_description=readme(),
      packages=find_packages('src'),
      package_dir={'': 'src'},
      install_requires=['networkx', 'six'],
      setup_requires=['pytest-runner', ],
      tests_require=['pytest', 'mock'],
      ext_modules=extensions,
      cmdclass={'build_ext': build_ext_new},
      zip_safe=False)
