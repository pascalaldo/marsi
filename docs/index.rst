=================
Welcome to marsi!
=================

|PyPI| |License| |Build Status| |Coverage Status|

**marsi** is an open-source software to created to identify non-GMO strain design targets.

There are two main experimental scenarios:

1. Adaptive Laboratory Evolution (ALE)
2. Classic Strain Improvement (CSI)

CSI
===

In this scenario we assume that metabolizing the target compound is going to kill the cells. Using chemical
mutagenesis, the surviving cells have found a way around the metabolism and are capable of resuming their activity
without the reactions related with that metabolite.

Here the search can be performed using existing methods (such as OptGene[1] or OptKnock[2]) that can predict knockout
targets. The targets will be then replaced and tested for the presence of an analog. We also implemented OptMet, a new
method that uses Heuristic Optimization to search for metabolite targets directly.


ALE
===

The ALE scenario assumes a long term exposition and adaptation of the cells to an analog metabolite. Here we account for
essential and non-essential metabolites. For essential metabolites the cells will produce more of the target metabolite
so it can compete with the analog for the enzymes. For non-essential metabolites we assume reduced activity/specificity
towards the target metabolite and the activity will be inhibited. To identify wich pathways should be inhibited we use
DifferentialFVA[3].



User's guide
============

.. toctree::
    :maxdepth: 1

    installation

.. toctree::

    examples/examples

Command Line Interface (CLI)
============================

.. toctree::
    :maxdepth: 1

    cli

Application Programming Interface (API)
=======================================

.. toctree::
    :maxdepth: 2

    API

References
==========
[1] Patil,K.R. et al. (2005) Evolutionary programming as a platform for in silico metabolic engineering. BMC Bioinformatics, 6, 308.

[2] Burgard,A.P. et al. (2003) Optknock: a bilevel programming framework for identifying gene knockout strategies for microbial strain optimization. Biotechnol. Bioeng., 84, 647–657.

[3] Cardoso,J.G.R. et al. (2017) Cameo : A Python Library for Computer Aided Metabolic Engineering and Optimization of Cell Factories. bioRxiv.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |PyPI| image:: https://img.shields.io/pypi/v/marsi.svg
   :target: https://pypi.python.org/pypi/marsi
.. |License| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: http://img.shields.io/badge/license-APACHE2-blue.svg
.. |Build Status| image:: https://travis-ci.org/biosustain/marsi.svg?branch=master
   :target: https://travis-ci.org/biosustain/marsi
.. |Coverage Status| image:: https://img.shields.io/codecov/c/github/biosustain/marsi/master.svg
   :target: https://codecov.io/gh/biosustain/marsi/branch/master
