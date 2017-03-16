from setuptools import setup


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
  version = '3.0.23',
  description = 'POVME (Pocket VOlume MEasurer) is a Python package for extracting actionable information from ensembles of protein structures for use in drug design.',
  author = 'Jeff Wagner',
  author_email = 'jwagnerjpl@gmail.com',
  url = 'https://github.com/POVME/POVME',
  download_url = 'https://github.com/POVME/POVME/tarball/2.0',
  setup_requires=['numpy','scipy'],
  #build_requires=['scipy','numpy'],
  #install_requires=[''],
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
