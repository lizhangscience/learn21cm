'''
Created on 09 May 2017
@author: Sambit Giri
Setup script
'''

from setuptools import setup, find_packages
setup(name='learn21cm',
      version='0.1',
      author='Sambit Giri',
      author_email='sambit.giri@astro.su.se',
      package_dir = {'learn21cm' : 'src'},
      packages=['learn21cm'],
      include_package_data=True,
)
