marsi - Metabolite Analogues for Strain Improvement
===================================================

marsi is an open-source software to created to identify non-GMO strain design
targets 

Dependencies
------------
* eigen3 >= 3.0
* OpenBabel >= 2.2.3
* RDKit >= 2016
* mongodb
* glpk
* cplex (optional)
* Cython >= 0.25
* numpy >= 1.11


Quick Start
-----------

1. Install the Dependencies
2. `pip install marsi` 
3. Run `marsi --help` from the command line

More details in [Quick Start](QUICK_START.md)

Initialization
--------------

*marsi* comes with a initialization command that will download all the necessary files
and build the database. You can start by running `marsi init --help`.

*marsi* will download the required files for you. Just run `marsi init download`. **Make sure you have an stable internet 
connection.** 


The PubChem API does not provide a method for query. To retrieve the necessary files,
go to [PubChem](https://pubchem.ncbi.nlm.nih.gov) and enter the following query:
 ```
    (antimetabolites) OR (analog) OR (analogue)     
 ```
 
And download the summary file. You should save the file in data folder at *marsi*'s working 
directory (`$HOME/.marsi/data` by default) with the name `pubchem_compound_analogs_antimetabolites.txt`.


Find Metabolite Targets
-----------------------


License
-------
Apache License Version 2.0


Databases and Software Tools
----------------------------

All tools and databases are listed [here](CITATIONS.md)