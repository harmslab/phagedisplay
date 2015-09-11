#!/usr/bin/env python3

# Try using setuptools first, if it's installed
from setuptools import setup, find_packages

    
# Need to add all dependencies to setup as we go!
setup(name='phagedisplay',
      version='0.1',
      description='Sequence Space in the Protein Lattice Model Landscape',
      author='Luke Wheeler, Zach Sailer, Andrea Loes, Dr. Michael J. Harms',
      url='https://github.com/harmslab/phagedisplay',
      packages=find_packages(),
#[
#            "phagedisplay",
#            "phagedisplay.processors",
 #           "phagedisplay.simulate",
#            "phagedisplay.cluster"
#            ],
      zip_safe=False,
      install_requires=["scipy","numpy"])

