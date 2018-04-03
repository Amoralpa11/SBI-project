from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))



setup(name='pmcb',
      version='1.0',
      description='Protein Macro-Complex Builder',
      url='https://github.com/Amoralpa11/SBI-project/tree/complex_breaker',
      author='Marc Elosua Bayés, Adrián Morales Pastor',
      author_email='elosua.marc@gmail.com, drnmoralespastor@gmail.com',
      packages=find_packages(),
      entry_points={'console_scripts': ['build_complex = build_complex:main']},
      install_requires=['biopython', 'matplotlib', 'pandas'])
