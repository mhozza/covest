covest
======

Tool that estimates coverage (and genome size) of dna sequence from
reads.

Installation
------------

``pip install covest``

For development:
~~~~~~~~~~~~~~~~

``pip install -e .`` from the project directory

Usage
-----

type ``covest --help`` for the usage.

Other included tools
--------------------

-  ``geset.py`` tool for estimation genome size from reads and known
   coverage
-  ``prepare_experient.py`` tool for experiment pipeline setup
-  ``experiment_table.py`` tool which collects data from experiment and
   create a nice table (html, tex, and csv formats are supported)
-  ``sam_to_fasta.py`` tool for converting sam file to fasta
-  Read simulator:

   -  ``generate_sequence.py`` random sequence generator
   -  ``read_simulator.py`` tool for generating random reads form the
      sequence

Copyright and citation
----------------------

Covest is licenced under `GNU
GPLv3 <http://www.gnu.org/licenses/gpl-3.0.en.html>`__ license.

Covest is research software, so you should cite us when you use it in scientific publications!
   Hozza, M., Vinař, T., & Brejová, B. (2015, September). How Big is that Genome? Estimating Genome Size and Coverage from k-mer Abundance Spectra. In String Processing and Information Retrieval (pp. 199-209). Springer International Publishing.
