# Copyright 2016 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.
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
from Cython.Build import cythonize
from setuptools import setup, find_packages, Extension

import versioneer

requirements = []
# ['pandas>=1.4.2',
#                 'numpy>=1.22.3',
#                 'bioservices>=1.8.4',
#                 'requests>=2.27.1',
#                 'bokeh>=2.4.2',
#                 'cameo>=0.13.6',
#                 'sqlalchemy>=1.4.35',
#                 'scikit-learn',
#                 'psycopg>=3.0.11',
#                 'bitarray>=2.4.1',
#                 'Cython>=0.29.28',
#                 'cement>=3.0.6',
#                 'pubchempy>=1.0.4 ',
#                 'cachetools>=5.0.0',
#                 'alembic>=1.7.7',
#                 'gnomic>=1.0.1',
#                 'markupsafe==2.0.1',
#                 'plotly',
#                 'six']

extra_requirements = {
    'docs': ['Sphinx>=1.3.5', 'numpydoc>=0.5'],
    'jupyter': ['jupyter>=1.0.0', 'ipywidgets>=4.1.1'],
    'test': ['pytest>=1.3.7', 'pytest-cov>=2.4', 'pytest-benchmark>=3.0'],
    '3d': ['imolecule>=0.1.13'],
    'opencl': ['pyopencl>=2016.1']
}

extra_requirements['all'] = sum([list(values) for values in extra_requirements.values()], [])


ext_modules = cythonize([Extension("marsi.chemistry.common_ext", 
                                   sources=["marsi/chemistry/common_ext.pyx"],
                                   include_dirs=[numpy.get_include()]), 
                         Extension("marsi.nearest_neighbors.model_ext",
                                   sources=["marsi/nearest_neighbors/model_ext.pyx"],
                                   include_dirs=[numpy.get_include()])
                        ])

include_dirs = [numpy.get_include()]


setup(
    name='marsi',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(exclude=("*.c",)),
    package_data={'marsi': ['alembic.ini', 'alembic/*',
                            'chemistry/common_ext.pxd', 'chemistry/common_ext.pyx',
                            'nearest_neighbors/model_ext.pyx']},
    install_requires=requirements,
    extras_require=extra_requirements,
    ext_modules=ext_modules,
    # scripts=['bin/marsi'],
    include_package_data=True,
    author='Joao Cardoso',
    author_email='joaca@biosustain.dtu.dk',
    description='marsi - Metabolite Analogues for Rational Strain Improvement',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics chemoinformatics',
    url='https://github.com/biosustain/marsi',
    long_description="marsi is an open-source software to created to identify non-GMO strain design targets",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

    ],
    entry_points={
          'console_scripts': ['marsi = marsi.cli.app:main'],
    },
    include_dirs=include_dirs
)
