============
Installation
============

Setting up a virtual environment first
======================================

We highly recommended installing marsi inside using (anaconda_).

There is no conda recipe for marsi, but the dependencies can be easily installed using ``conda``.
*Marsi* needs to install it using ``pip``. Do the following to create a virtual environment and get
 some of the heavier dependencies out of the way.

.. code-block:: bash

    $ conda create -y -n marsi python=3.4 cython

Alternatively you can use ``virtualenv``. virtualenvwrapper_ tremendously simplifies using virtualenv_ and can easily
be installed using virtualenv-burrito_. Once you installed virtualenv_ and virtualenvwrapper_, run

.. code-block:: bash

    $ mkvirtualenv marsi  # or whatever you'd like to call your virtual environment
    $ workon marsi

and then continue with the installation instructions described below.


Dependencies required before installation
=========================================

Using conda
-----------

All dependencies will by installed by running the following commands

.. code-block:: bash

    $ conda install -c openbabel openbabel
    $ conda install -c rdkit rdkit


Alternativly you can install the dependencies one by one.

Cython
------
To install Cython follow the cython-installation_ documentation provided on their website.

Numpy
-----

Numpy can be installed using ``conda``:

.. code-block:: bash

    $ conda install numpy

Or can also be installed using pip:

.. code-block:: bash

    $ pip install numpy

Eigen3
------

To install Eigen3 on Ubuntu/Debian run:

.. code-block:: bash

    $ sudo apt-get install libeigen3-dev

To install Eigen3 on MacOS X use homebrew_:

.. code-block:: bash

    $ brew install eigen3


OpenBabel
---------

To install OpenBabel_ on Ubuntu/Debian run:

.. code-block:: bash

    $ sudo apt-get install openbabel libopenbabel-dev

On MacOS X it can be installed using homebrew_:

.. code-block:: bash

    $ brew install open-babel


RDKit
-----

RDKit_ can be installed using the rdkit-documentation_.

For MacOS X we recommend using homebrew_:

.. code-block:: bash

    $ brew install rdkit


GLPK
----

Using constraint-based methods with marsi requires glpk_ to be installed.
In order to generate python bindings, swig_ is also required.
On Ubuntu/Debing  we recommend using:

.. code-block:: bash

    $ sudo apt-get install libglpk-dev glpk-utils swig

On MacOS X it can be installed using homebrew_.

.. code-block:: bash

    $ brew install swig
    $ brew install glpk


Optional
--------

CPLEX
-----

CPLEX is a commercial solver with great performance. It solves linear, mixed-integer and quadratic problems.
Install CPLEX as described in cplex-install_ and then install the python bindings (cplex-python_).

Installation
============

**marsi** can be installed using ``pip``.

.. code-block:: bash

    $ pip install marsi


Soft dependencies
=================

The following soft dependencies can be installed all at once using ``pip install marsi[all]`` or individually
by specifying individual categories of dependencies (for example ``pip install marsi[jupyter,3d, ...]``).
The following categories are available::

    'docs': ['Sphinx>=1.3.5', 'numpydoc>=0.5'],
    'jupyter': ['jupyter>=1.0.0', 'ipywidgets>=4.1.1'],
    'test': ['nose>=1.3.7', 'rednose>=0.4.3', 'coverage>=4.0.3'],
    '3d': ['imolecule>=0.1.13'],
    'opencl': ['pyopencl>=2016.1']



.. _anaconda: https://anaconda.org
.. _homebrew: http://brew.sh/
.. _RDKit: http://www.rdkit.org
.. _OpenBabel: http://openbabel.org
.. _rdkit-documentation: http://www.rdkit.org/docs/Install.html
.. _glpk: https://www.gnu.org/software/glpk/
.. _cplex: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
.. _virtualenv-burrito: https://github.com/brainsik/virtualenv-burrito
.. _virtualenv: https://pypi.python.org/pypi/virtualenv
.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrapper
.. _cython-installation: http://cython.readthedocs.io/en/latest/src/quickstart/install.html
.. _sphinx: https://pypi.python.org/pypi/sphinx
.. _swig: http://www.swig.org
.. _numpydoc: https://pypi.python.org/pypi/numpydoc
.. _cplex-install: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.7.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html
.. _cplex-python: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html

