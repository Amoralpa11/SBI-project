from distutils.core import setup

setup(name='pmcb',
      version='1.0',
      description='Protein Macro-Complex Builder',
      url='https://github.com/Amoralpa11/SBI-project/tree/complex_breaker',
      author='Marc Elosua Bayés, Adrián Morales Pastor',
      author_email='elosua.marc@gmail.com, drnmoralespastor@gmail.com',
      packages=['pmcb'],
      scripts=['build_complex.py'],
      requires=['Bio', 're', 'os', 'copy', 'modeller', 'matplotlib', 'pandas'])
