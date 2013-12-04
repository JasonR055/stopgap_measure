#!/usr/bin/env python

from distutils.core import setup
import sys

try:
    import pandas
    #print "Found pandas version " + pandas.__version__ + "."
except ImportError:
    sys.exit("install requires: 'pandas >=0.8.0'\n \
              http://pandas.pydata.org/pandas-docs/stable/install.html")

try:
    import pysam
    #print "Found pysam version " + pysam.__version__ + "."
except ImportError:
    sys.exit("install requires: 'pysam'\n \
              http://code.google.com/p/pysam/")

classifiers = """
Development Status :: 2 - Alpha
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""

setup(
    name='stopgap_measure',
    version='0.1',
    description='Methylation estimation of bisulphite-treated DNA',
    author='Jason Ross',
    author_email='jason.ross@csiro.au',
    url='http://www.csiro.au/fns',
    license='CSIRO license <http://www.ict.csiro.au/downloads.php>',
    packages=['measure'],
    platforms='ALL',
    scripts = ['run_measure.py'],
    requires=['pandas (>=0.8.0)', 'pysam']
    )

