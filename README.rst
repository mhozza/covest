CovEst
======

Tool that estimates coverage (and genome size) of dna sequence from
reads.

.. image:: https://badge.fury.io/py/covest.svg
    :target: https://badge.fury.io/py/covest
.. image:: https://travis-ci.org/mhozza/covest.svg?branch=master
    :target: https://travis-ci.org/mhozza/covest

Requirements
------------
- python 3.4+
- python3-dev
- gcc

Installation
------------
We suggest to install *CovEst* in python3 virtual environment.

``pip install covest``

For development:
~~~~~~~~~~~~~~~~

``pip install -e .`` from the project directory

Usage
-----

type ``covest --help`` for the usage.

Basic Usage:
~~~~~~~~~~~~
``covest histogram -m model -k K -r read_length``

-  You can specify the read file using ``-s reads.fa`` parameter for more precise genome size computation.
-  default *K* is 21
-  default *read length* is 100
-  currently, the supported models are:

   -  basic: for simple genomes without repeats
   -  repeat: for genomes with repetitive sequences

Input Histogram Specification:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The input histogram can be generated from the read data using `jellyfish <http://www.cbcb.umd.edu/software/jellyfish/>`__.

-  ``jellyfish count -m K -C reads.fa -o table.jf``
-  ``jellyfish histo table.jf -o reads.hist``

The format of the histogram is just list of lines. Each lines contains an index and value separated by space.

Output Specification:
~~~~~~~~~~~~~~~~~~~~~
CovEst outputs it's results in simple subset of YAML format for best human readability and possibility of machine processing.

The output are lines containing ``key: value``. The most important keys are ``coverage`` and ``genome_size`` (or ``genome_size_reads`` if read file was specified).

Other included tools
--------------------

-  ``geset.py`` tool for estimation genome size from reads and known
   coverage
-  ``kmer_hist.py`` custom khmer histogram computation, it is much slower than other tools, so use it only if you have no other option.
-  ``read_sampler.py`` script for subsampling reads, useful if you have very high coverage data and want to make it smaller.
-  ``fasta_length.py`` get total length of all sequences in fasta file.
-  Read simulator:

   -  ``generate_sequence.py`` random sequence generator
   -  ``read_simulator.py`` tool for generating random reads form the
      sequence

Copyright and citation
----------------------

CovEst is licenced under `GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.en.html>`__ license.

CovEst is research software, so you should cite us when you use it in scientific publications!
   Hozza, M., Vinař, T., & Brejová, B. (2015, September). How Big is that Genome? Estimating Genome Size and Coverage from k-mer Abundance Spectra. In String Processing and Information Retrieval (pp. 199-209). Springer International Publishing.
