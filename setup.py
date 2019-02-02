try:
    from setuptools import setup
except:
    from distutils.core import setup


import re, ast
_version_re = re.compile(r'__version__\s+=\s+(.*)')
with open('hpic/__init__.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

long_description = '''Pure Ion Chromatogram Extraction via Hierarchical Density Clustering for LC-MS.'''

setup(name='hpic',
      version=version,
      description='Pure Ion Chromatogram Extraction via Hierarchical Density Clustering for LC-MS.',
      long_description=long_description,
      author='Zhimin Zhang, Rong Wang',
      author_email='zhangzhimin.csu@gmail.com',
      license='MIT',
      url='https://github.com/zmzhang/HPIC',
      download_url='https://github.com/zmzhang/HPIC',
      packages=['hpic', 'hpic.fileio', 'hpic.mspd','hpic.hpic'],
      platforms=['Platform-Independent'],
      install_requires=[
	  "numpy>=1.15.0",
	  "scipy>=1.2.0",
	  "pyopenms>=2.4.0",
	  "hdbscan>=0.8.0"],
      classifiers=[ 
	    "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Metabolomics",
        "Topic :: Scientific/Engineering :: Chemistry",
      ]
)