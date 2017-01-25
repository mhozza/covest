#!/usr/bin/env python
# -*- coding: utf-8 -*-
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import sys

try:
    with open('README.rst') as readme_file:
        readme = readme_file.read()
except FileNotFoundError:
    readme = ''

requirements = [
    'biopython',
    'matplotlib',
    'numpy',
    'pystache',
    'scipy>=0.18',
    'pyfasta',
    'first',
    'pyyaml',
]

if sys.version_info < (3, 5):
    requirements.append('pyparsing!=2.1.2')

test_requirements = [
]

covest_poisson = Extension(
    'covest_poisson',
    libraries=['m'],
    sources=['c_src/covest_poissonmodule.c'],
    extra_compile_args=['-std=c99'],
)

setup(
    name='covest',
    version='0.5.6',
    description="Covest estimates the coverage and genome size, "
    "just from k-mer abundance histogram computed from DNA sequences reads.",
    long_description=readme,
    author="Michal Hozza",
    author_email='mhozza@gmail.com',
    url='https://github.com/mhozza/covest',
    packages=[
        'covest',
    ],
    package_dir={
        'covest': 'covest'
    },
    include_package_data=True,
    ext_modules=[covest_poisson],
    install_requires=requirements,
    entry_points={
        'console_scripts': ['covest=covest.covest:run'],
    },
    scripts=[
        'bin/gsest.py',
        'bin/reads_size.py',
        'bin/fasta_length.py',
        'bin/kmer_hist.py',
        'bin/read_sampler.py',
    ],
    license='GNU GPLv3',
    zip_safe=False,
    keywords='covest',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
