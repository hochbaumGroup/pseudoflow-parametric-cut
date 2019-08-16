# Copyright © 2001. The Regents of the University of California (Regents). All Rights
# Reserved.
#
# Permission to use, copy, modify, and distribute this software and its documentation
# for educational, research, and not-for-profit purposes, without fee and without a
# signed licensing agreement, is hereby granted, provided that the above copyright
# notice, this paragraph and the following two paragraphs appear in all copies,
# modifications, and distributions. Contact The Office of Technology Licensing, UC
# Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201,
# for commercial licensing opportunities. Created by Created by Bala Chandran, Quico
# Spaen, and Dorit S. Hochbaum, Department of Industrial Engineering and Operations
# Research, University of California, Berkeley.
#
# IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
# INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE
# OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
# SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
# IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
# ENHANCEMENTS, OR MODIFICATIONS.

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
            return ext_name.split(".")[-1] + ".so"
        return build_ext.get_ext_filename(self, ext_name)


class CTypes(Extension):
    pass


with open("README.md") as f:
    readme = f.read()


extensions = [
    CTypes(
        "pseudoflow.libhpf",
        ["src/pseudoflow/core/libhpf.c"],
        depends=["src/pseudoflow/core/libhpf.h"],
        export_symbols=["hpf_solve", "libfree"],
        # include_dirs=["pseudoflow/core"],
        language="c99",
        extra_compile_args=["-std=c99", "-O3"],
    )
]


setup(
    name="pseudoflow",
    version="2019.8.1",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
    ],
    description="Pseudoflow algorithm for the parametric minimum cut problem.",
    keywords=["minimum cut", "network flow", "parametric"],
    url="https://github.com/quic0/pseudoflow",
    author="Quico Spaen",
    author_email="qspaen@berkeley.edu",
    long_description_content_type="text/markdown",
    license="Non-commercial license. Not an open-source license.",
    long_description=readme,
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=["six"],
    ext_modules=extensions,
    cmdclass={"build_ext": build_ext_new},
    zip_safe=False,
)
