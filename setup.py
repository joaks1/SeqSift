#! /usr/bin/env python

from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='SeqSift',
      version=version,
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
      include_package_data=True,
      zip_safe=True,
      test_suite="seqsift.test",
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
