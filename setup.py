#!/usr/bin/env python

import os
import sys
import glob

# Try and import pip. We'll stop if it is not present
try:
    import pip
except ImportError:
    print "Installation of Contiguity requires pip. Please install it! See -"
    print "http://pip.readthedocs.org/en/latest/installing.html"
    sys.exit(1)

from setuptools import setup

__title__ = 'Contiguity'
__version__ = '1.0.2'
__description__ = "Tool for visualising assemblies"
__author__ = 'Mitchell Sullivan'
__license__ = 'GPLv3'
__author_email__ = "mjsull@gmail.com"
__url__ = 'https://github.com/BeatsonLab-MicrobialGenomics/Contiguity'


# Helper functions
if sys.argv[-1] == 'publish':
    print "Please use twine or do_release.sh"
    sys.exit()

if sys.argv[-1] == 'clean':
    os.system('rm -rf Contiguity.egg-info build dist')
    sys.exit()

if sys.argv[-1] == 'docs':
    os.system('cd docs && make html')
    sys.exit()


packages = [__title__, ]

requires = []
with open('requirements.txt') as fin:
    lines = fin.readlines()
    for line in lines:
        requires.append(line.strip())

# Build lists to package the docs
html, sources, static = [], [], []
html_f = glob.glob('docs/_build/html/*')
accessory = glob.glob('docs/_build/html/*/*')
for f in html_f:
    if os.path.isfile(f):
        html.append(f)
for f in accessory:
    if f.find("_static") != -1:
        if os.path.isfile(f):
            static.append(f)
    elif f.find("_sources"):
        if os.path.isfile(f):
            sources.append(f)

setup(
    name=__title__,
    version=__version__,
    description=__description__,
    long_description=open('README.rst').read(),
    author=__author__,
    author_email=__author_email__,
    url=__url__,
    packages=packages,
    test_suite="tests",
    package_dir={__title__.lower(): __title__},
    scripts=[__title__.lower()+'/'+__title__],
    package_data={},
    data_files=[('', ['LICENSE', 'requirements.txt', 'README.rst']),
                ('docs', html), ('docs/_static', static),
                ('docs/_sources', sources)],
    include_package_data=True,
    install_requires=requires,
    license=__license__,
    zip_safe=False,
    classifiers=('Development Status :: 4 - Beta',
                 'Environment :: X11 Applications',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved',
                 'Natural Language :: English',
                 'Operating System :: POSIX :: Linux',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 2 :: Only',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Visualization',),
)
