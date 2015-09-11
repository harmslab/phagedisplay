#!/usr/bin/env python3

import sys

# Try using setuptools first, if it's installed
from setuptools import setup, find_packages


if not sys.version_info[0] == 3:
    sys.exit("Cannot set up.  Python 3 is required.")
    
# Need to add all dependencies to setup as we go!
setup(name='phagedisplay',
      version='0.1.1',
      description='Sequence Space in the Protein Lattice Model Landscape',
      author='Luke Wheeler, Zach Sailer, Andrea Loes, Dr. Michael J. Harms',
      url='https://github.com/harmslab/phagedisplay',
      packages=find_packages(),
      zip_safe=False,
      install_requires=["scipy","numpy"])

