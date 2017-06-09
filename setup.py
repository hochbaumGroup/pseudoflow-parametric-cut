from setuptools import setup, Extension
import unittest


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
    Extension(
        "pseudoflow.libhpf",
        ["pseudoflow/core/libhpf.c"],
        depends=["pseudoflow/core/libhpf.h"],
        include_dirs=["pseudoflow/core"],
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
      keywords='minimum cut network flow combinatorial algorithms parametric',
      url='https://github.com/quic0/pseudoflow',
      author='Quico Spaen',
      author_email='qspaen@berkeley.edu',
      license=license(),
      long_description=readme(),
      packages=['pseudoflow', 'pseudoflow/python'],
      install_requires=['networkx', 'six'],
      test_suite='setup.test_suite',
      ext_modules=extensions,
      zip_safe=False)
