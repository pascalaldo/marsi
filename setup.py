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

from setuptools import setup, find_packages

requirements = ['pandas>=0.18.1',
                'numpy>=1.11.1',
                'bioservices>=1.4.14',
                'requests>=2.11.1',
                'bokeh>=0.12',
                'cameo>=0.7']

extra_requirements = {}

setup(
    name='marsi',
    version="0.0.1a",
    packages=find_packages(),
    install_requires=requirements,
    extras_require=extra_requirements,
    include_package_data=True,
    author='Joao Cardoso',
    author_email='joaca@biosustain.dtu.dk',
    description='marsi - Metabolite Analogues for Rational Strain Improvement',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics',
    url='http://cameo.bio',
    long_description="",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.5',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License'
    ],
)
