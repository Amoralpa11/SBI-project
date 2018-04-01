from distutils.core import setup

setup(name='pmcb',
      version='1.0',
      description='Protein Macro-Complex Builder',
      author='Marc Elosua Bayés, Adrián Morales Pastor',
      author_email='elosua.marc@gmail.com, drnmoralespastor@gmail.com',
      py_modules=["modules/Complex_breaker", "modules/Complex_id", "modules/macrocomplex_builder", 'modules/modeller_optimization',
                  'modules/pdb_files_comparison'],
      scripts=['build_complex.py'])
