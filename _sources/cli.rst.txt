CLI Docs
========

*marsi* provides a command line interface (CLI) that allows running *marsi* from the command line.

The CLI is split in three parts.

`marsi db`

`marsi optimize`

`marsi chem`

marsi db
--------

Contains the commands to build the chemical analogues database.

To get started type:

.. code-block:: bash

    $ marsi db init


This will take a few hours to download the necessary files.

marsi optimize
--------------

The commands under this subcommand are used to perform constraint-based modeling.

There are two main functions:

.. code-block:: bash

    $ marsi optimize mutagenesis

will calculate knockout-based designs. And

.. code-block:: bash

    $ marsi optimize ale

will calculate knockout based and over-, down-regulated designs.

marsi chem
----------

The `chem` subcommand can be used to search the chemical analogues database.

.. code-block:: bash

    $ marsi chem find-analogs --inchi=<inchi-str>
