from distutils.core import setup
from Cython.Build import cythonize
import numpy as np


setup(
    name = 'DistMatrix',
    ext_modules = cythonize('DistMatrix.pyx'),
    include_dirs = [np.get_include()]
)