#!/usr/bin/env python

from setuptools import setup
import glob

from eccFinder_lib.utilities import get_eccFinder_version

scripts = glob.glob("*.p*")

setup(
    name='ecc_finder',
    version='v1.0.0',
    description='A tool for detecting extrachromosomal circular DNA (eccDNA) from sequencing data.',
    author='Panpan Zhang',
    author_email='njaupanpan@gmail.com',
    packages=['eccFinder_lib'],
    package_dir={'eccFinder_lib': 'eccFinder_lib/'},
    install_requires=[
              'pysam',
              'numpy',
              'pandas',
              'matplotlib',
              'pybedtools'
          ],
    scripts=scripts,
    zip_safe=True
