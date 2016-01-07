from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_module = Extension(
    "DistMatrix",
    ["DistMatrix.pyx"],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp']
)

setup(
    name = 'Hello world app',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_module],
    include_dirs = [np.get_include()]
)