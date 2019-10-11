from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module1 =  Extension('Bpowspec',
                     sources = ['code/Bpowspec.c','code/Bpowspec_mod.c'],
                     include_dirs = [numpy_inc,'code'],
                     libraries=['gsl','gslcblas','fftw3'],
                     library_dirs = ["lib"],
                     extra_compile_args=['-fPIC','-Wall','-g'])

setup (name = 'Bpowspec',
       version = '0.1',
       url='https://github.com/huffenberger-cosmology/magnetic-field-power-spectrum',
       description = '',
       ext_modules = [module1],
       )

