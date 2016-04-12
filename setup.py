#!/usr/bin/env python3

import sys

# Try using setuptools first, if it's installed
from setuptools import setup, find_packages


if not sys.version_info[0] == 3:
    sys.exit("Cannot set up.  Python 3 is required.")
    
# Need to add all dependencies to setup as we go!
setup(name='phagedisplay',
      version='0.1.1',
      description='Pipeline for quantitative analysis of phage display experiments',
      author='Luke Wheeler, Hiranmyai Duvvuri, Michael J. Harms',
      url='https://github.com/harmslab/phagedisplay',
      packages=find_packages(),
      zip_safe=False,
      install_requires=["scipy","numpy","pandas","sklearn"])

