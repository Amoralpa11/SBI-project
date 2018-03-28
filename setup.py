from distutils.core import setup

setup(name='pmcb',
      version='1.0',
      description='Protein Macro-Complex Builder',
      author='Marc Elosua Vallès, Adrián Morales Pastor',
      author_email='elosua.marc@gmail.com, drnmoralespastor@gmail.com',
      py_modules=["Complex_breaker", "Complex_id", "macrocomplex_builder", 'modeller_optimization', 'pdb_files_comparison'])