# Copyright 2016 Chr. Hansen A/S.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import absolute_import, print_function

import numpy
from setuptools import setup, find_packages
from Cython.Build import cythonize

requirements = ['pandas>=0.18.1',
                'numpy>=1.11.1',
                'bioservices>=1.4.14',
                'requests>=2.11.1',
                'bokeh>=0.12',
                'cameo>=0.7.1',
                'openbabel>=1.8.4',
                'mongoengine>=0.10.6',
                'bitarray>=0.8.1',
                'Cython>=0.24',
                'cement>2.10']

extra_requirements = {
    'docs': ['Sphinx>=1.3.5', 'numpydoc>=0.5'],
    'jupyter': ['jupyter>=1.0.0', 'ipywidgets>=4.1.1'],
    'test': ['nose>=1.3.7', 'rednose>=0.4.3', 'coverage>=4.0.3'],
    '3d': ['imolecule>=0.1.13'],
    'opencl': ['pyopencl>=2016.1']
}

extra_requirements['all'] = sum([list(values) for values in extra_requirements.values()], [])


ext_modules = cythonize(["marsi/processing/chemistry/common_ext.pyx", "marsi/mining/_nearest_neighbors_ext.pyx"],
                        include_dirs=[numpy.get_include()])

include_dirs = [numpy.get_include()]


setup(
    name='marsi',
    version="0.0.1a0",
    packages=find_packages(exclude=("*.pyx")),
    install_requires=requirements,
    extras_require=extra_requirements,
    ext_modules=ext_modules,
    scripts=['bin/marsi'],
    include_package_data=True,
    author='Joao Cardoso',
    author_email='joaca@biosustain.dtu.dk',
    description='marsi - Metabolite Analogues for Rational Strain Improvement',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics chemoinformatics',
    url='http:/nourl.com',
    long_description="",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.4',
        'License :: OSI Approved :: Apache Software License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

    ],
    include_dirs=include_dirs
)
