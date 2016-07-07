from distutils.core import setup

setup(
  name = 'POVME',
  packages = ['POVME'],
  version = '3.0',
  description = 'POVME (Pocket VOlume MEasurer) is a Python package for extracting actionable information from ensembles of protein structures for use in drug design.',
  author = 'Jeff Wagner',
  author_email = 'jwagnerjpl@gmail.com',
  url = 'https://github.com/POVME/POVME',
  download_url = 'https://github.com/POVME/POVME/tarball/2.0',
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
