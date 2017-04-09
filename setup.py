from setuptools import setup

setup(name='quantitative_mri',
      version='0.1',
      description='A collection of functions to analyze CEST MRI / NMR data',
      url='http://github.com/JCardenasRdz/cest_mri',
      author='Julio Cardenas-Rodriguez',
      author_email='jdatscientist@gmail.com',
      license='Apache',
      packages=['quantitative_mri'],
      install_requires=['numpy', 'scipy'],
      zip_safe=False)
