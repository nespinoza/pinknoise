import re
from setuptools import setup

VERSIONFILE='pinknoise/_version.py'
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(name='pinknoise',
      version=verstr,
      description='pinknoise: a library to perform inference and simulations of all flavors of pink noise',
      url='http://github.com/nespinoza/pinknoise',
      author='Nestor Espinoza',
      author_email='nespinoza@stsci.edu',
      license='MIT',
      packages=['pinknoise'],
      install_requires=['numpy','scipy'],
      python_requires='>=3.0',
      extras_requires={
            'seaborn':['seaborn'],
            'matplotlib':['matplotlib'],},
      zip_safe=False)
