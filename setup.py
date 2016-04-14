#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'biopython',
    'matplotlib',
    'numpy',
    'pystache',
    'scipy<0.16',
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='covest',
    version='0.1.0',
    description="Covest estimates the coverage and genom size, "
    "just from k-mer abundance histogram computed from DNA sequences reads.",
    long_description=readme + '\n\n' + history,
    author="Michal Hozza",
    author_email='mhozza@gmail.com',
    url='https://github.com/mhozza/covest',
    packages=[
        'covest',
    ],
    package_dir={'covest':
                 'covest'},
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': ['covest=covest.covest:run'],
    },
    license='ISCL',
    zip_safe=False,
    keywords='covest',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
