'''
A setuptools based setup module to install willowtree on your personal
computer.

How does it work?
-------------------------------------------------------------------------------
On Terminal (macOS) or Command Prompt (Windows), navigate through folder
'willow-tree', then run ($ is the Terminal prompt):

$ python3 setup.py install

This command will generate three folders in the current directory:

build
dist
willowtree.egg-info

and one .egg file in the main package folder, for example, if third-party
extension Anaconda is installed, in 'anaconda/lib/python3.X/site-packages':

willowtree-X.X-py3.X.egg

where X.X refers to the most recent installed version of 'willowtree', e.g.
0.9, and py3.X to the latest installed Python version, e.g. 3.6.

It is now possible to import willowtree in a Python environment.

References
-------------------------------------------------------------------------------
https://pythonhosted.org/an_example_pypi_project/setuptools.html
https://github.com/pypa/sampleproject
'''

from setuptools import setup

# Look for module __version__.py in willowtree and import variable __version__
from willowtree.__version__ import __version__

# Assign the value of __version__ to variable version, then use it in setup()
version = __version__

setup(
    name='willowtree',

    version=version,

    description='''Robust and flexible Python implementation of the willow
    tree lattice for derivatives pricing.''',

    url='https://github.com/federicomariamassari/willow-tree',

    author='Federico Maria Massari',
    author_email='federico.massari@bocconialumni.it',

    license='MIT',

    classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Financial and Insurance Industry',
    'Topic :: Scientific/Engineering :: Mathematics',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'
    ],

    keywords=[
    'willow-tree',
    'derivatives-pricing',
    'standard-brownian-motion'
    ],

    # The packages to install. There is only one folder, 'willowtree'
    packages=[
    'willowtree'
    ],

    # Dependencies, the auxiliary libraries necessary to run 'willowtree'
    install_requires=[
    'numpy >= 1.13',
    'scipy >= 0.19',
    'matplotlib >= 2.0',
    'seaborn >= 0.8'
    ],
)
