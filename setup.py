import re
import numpy
#from distutils.core import setup, Extension
from setuptools import setup, Extension

VERSIONFILE='src/pinknoise/_version.py'
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

module = Extension('FWT', sources = ['src/pinknoise/FWT.c'], libraries=['m'], include_dirs=[numpy.get_include(),'/usr/local/include'])
setup(name='pinknoise',
      version=verstr,
      description='pinknoise: a library to perform inference and simulations of all flavors of pink noise',
      url='http://github.com/nespinoza/pinknoise',
      author='Nestor Espinoza',
      author_email='nespinoza@stsci.edu',
      license='MIT',
      packages=['pinknoise'],
      package_dir={'pinknoise': 'src/pinknoise'},
      install_requires=['numpy','scipy', 'stochastic'],
      python_requires='>=3.0',
      ext_modules = [module],
      zip_safe=False)
