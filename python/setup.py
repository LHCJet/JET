from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
extensions = [
    Extension("jetjet", ["jetjet.pyx"],
        include_dirs = ["../src/"],
        libraries = ["JETJet"],
        library_dirs = ['../build/'],
        language="c++",
        extra_compile_args=["-std=c++11"])
]

setup(
    name = "JET algorithm",
    ext_modules = cythonize(extensions),
)
