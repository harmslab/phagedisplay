from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

extensions = Extension(
    "DistMatrix",
    ["DistMatrix.pyx"]
)

setup(
    name = 'DistMatrix',
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(extensions),
    include_dirs = [np.get_include()]
)