from setuptools import setup, find_packages


setup(
  name = 'POVME',
  packages = find_packages(),
  scripts = ['POVME/packages/clustering/binding_site_overlap.py',
             'POVME/packages/clustering/cluster.py',
             'POVME/POVME3.py',
             'POVME/packages/clustering/pocketPointsPca.py'], 
  version = '3.0.34',
  description = 'POVME (Pocket VOlume MEasurer) is a Python package for extracting actionable information from ensembles of protein structures for use in drug design.',
  author = 'Jeff Wagner',
  author_email = 'j5wagner@ucsd.edu',
  url = 'https://github.com/POVME/POVME',
  download_url = 'https://github.com/POVME/POVME/tarball/2.0',
  install_requires=['scipy','numpy','matplotlib','networkx'],
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
