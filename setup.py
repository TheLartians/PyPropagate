from setuptools import setup, Extension, find_packages
from glob import glob
import numpy

setup(
    name='pypropagate',
    version='1.0.1',
    description='A python based paraxial wave propagation framework',
    
    author='Lars Melchior',
    author_email='lars.melchior@gmail.com',

    packages=find_packages(exclude=['tests*']),
    
    install_require={
        'matplotlib',
        'numpy',
        'expresso[pyplot]'
    },

    zip_safe=False,

    classifiers=[
        'Programming Language :: Python :: 2.7'
    ],

    ext_modules=[
        Extension('_pypropagate',
                  sources = ['source/crank_nicolson.cpp','source/finite_difference.cpp','source/python.cpp'],
                  include_dirs=['libs',numpy.get_include()], 
                  libraries=['boost_python'], 
                  library_dirs=['/'],
                  extra_compile_args=['-g','-std=c++11','-Wno-unknown-pragmas','-O3']
                  ),
        ]
)
