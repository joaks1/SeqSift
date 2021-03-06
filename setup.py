#! /usr/bin/env python

from setuptools import setup, find_packages
import sys, os
import glob

from seqsift import __version__, __project__

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
SCRIPTS = glob.glob(os.path.join(BASE_DIR, 'scripts', '*.py'))

setup(name=__project__,
      version=__version__,
      description="A python package for vetting DNA sequence data",
      long_description="""\
A python package for vetting DNA sequence data""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='DNA bioinformatics',
      author='Jamie Oaks',
      author_email='joaks1@gmail.com',
      url='',
      license='GPL',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=SCRIPTS,
      include_package_data=True,
      zip_safe=True,
      test_suite="seqsift.test",
      install_requires=[
          'biopython'
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
