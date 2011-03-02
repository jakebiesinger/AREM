#!/usr/bin/env python
# Time-stamp: <2011-03-01 18:21:42 Jake Biesinger>

"""Description:

Copyright (c) 2008,2009,2010 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status: beta
@version: $Revision$
@originalauthor:  Yong Zhang, Tao Liu
@originalcontact: taoliu@jimmy.harvard.edu

Modifications to probabilistically align reads to regions with highest
enrichment performed by Jacob Biesinger. Repackaged as "AREM" in accordance
with copyright restrictions.

@author: Biesinger, W Jacob B
@contact: jake.biesinger@gmail.com

Changes to this file since original release of MACS 1.4 (summer wishes):
  December/January 2011
    * Updated names (AREM, not MACS14)
"""

import os
import sys
try:
    from setuptools import setup, Extension
except:
    from distutils.core import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def main():
    if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
        sys.exit(1)

    setup(name="AREM",
          version="1.0.1",
          description="Aligning Reads by Expectation-Maximization.\nBased on MACS (Model Based Analysis for ChIP-Seq data)",
          author='Jake Biesinger; Daniel Newkirk; Alvin Chon; Yong Zhang; Tao (Foo) Liu',
          author_email='jake.biesinger@gmail.com; dnewkirk@uci.edu; achon@uci.edu; zy@jimmy.harvard.edu; taoliu@jimmy.harvard.edu',
          url='http://cbcl.ics.uci.edu/AREM',
          long_description=read('README'),
          package_dir={'AREM' : 'AREM'},
          packages=['AREM', 'AREM.IO'],
          scripts=['bin/arem','bin/elandmulti2bed','bin/elandresult2bed','bin/elandexport2bed',
                   'bin/sam2bed', 'bin/wignorm'],
          license = "Creative",
          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: Artistic License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Programming Language :: Python :: 2.6',
              'Programming Language :: Python :: 2.7',
              'Topic :: Scientific/Engineering :: Bio-Informatics'
              ],
          )

if __name__ == '__main__':
    main()
