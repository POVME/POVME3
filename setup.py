from setuptools import setup




#!/usr/bin/env python
# Install dependencies listed in pyproject.toml, then
# continue with regularly scheduled setup.py.
# From https://bitbucket.org/dholth/setup-requires

import sys, subprocess, pkg_resources

#sys.path[0:0] = ['setup-requires']
#pkg_resources.working_set.add_entry('setup-requires')

def missing_requirements(specifiers):
    for specifier in specifiers:
        try:
            pkg_resources.require(specifier)
        except pkg_resources.DistributionNotFound:
            yield specifier

def install_requirements(specifiers):
    to_install = list(specifiers)
    if to_install:
        subprocess.call([sys.executable, "-m", "pip", "install", 
            ] + to_install)
        #subprocess.call([sys.executable, "-m", "pip", "install", 
        #    "-t", "setup-requires"] + to_install)
        
install_requirements(missing_requirements(['pytoml']))

import pytoml

try:
    with open('pyproject.toml') as f:
        pyproject = pytoml.load(f)
except IOError:
    pass
else:
    requires = pyproject.get('build-system', {}).get('requires')
    install_requirements(missing_requirements(requires))

    '''
### Place normal setup.py contents below ###

from setuptools import setup


setup(name="example-package",
    version = "0.0.1",
    py_modules = [ 'example_package' ],
    install_requires = [ ],
    description = "An example of a package with setup requirements.",
    license = "MIT",
    author = "Emilio Example",
    author_email = "emilio@example.org",
    url="https://bitbucket.org/dholth/setup-requires")

# Alternatively, run real-setup.py from a separate file
# exec(compile(open("real-setup.py").read().replace('\\r\\n', '\\n'),
#     __file__,
#     'exec'))

'''









setup(
  name = 'POVME',
  packages = ['POVME',
              'POVME.packages',
              'POVME.packages.pymolecule',
              'POVME.packages.binana',
              'POVME.packages.clustering'],
  scripts = ['POVME/packages/clustering/binding_site_overlap.py',
             'POVME/packages/clustering/cluster.py',
             'POVME/POVME3.py',
             'POVME/packages/clustering/pocketPointsPca.py'],
  package_dir={'POVME': 'POVME'},
  package_data={'POVME': ['examples/*/*.*','examples/*/*/*.*'#,'POVME/examples/*/*/*/*'
                          ]},        
  version = '3.0.30',
  description = 'POVME (Pocket VOlume MEasurer) is a Python package for extracting actionable information from ensembles of protein structures for use in drug design.',
  author = 'Jeff Wagner',
  author_email = 'jwagnerjpl@gmail.com',
  url = 'https://github.com/POVME/POVME',
  download_url = 'https://github.com/POVME/POVME/tarball/2.0',
  #setup_requires=['scipy','numpy'],
  #build_requires=['scipy','numpy'],
  install_requires=['scipy','numpy'],
  keywords = ['protein folding', 'protein', 'pocket volume measurer', 'computational protein folding', 'protein structures', 'drug design'],
  classifiers = [
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Programming Language :: Python :: 2.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Chemistry"
  ]
)
