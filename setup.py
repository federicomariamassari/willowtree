'''
A setuptools bases setup module for the willow tree.

Also see: https://github.com/pypa/sampleproject
'''

from setuptools import setup

setup(
    name='willow-tree',

    version='0.9',

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

    packages=[
    'python-modules'
    ],

    install_requires=[
    'numpy >= 1.13',
    'scipy >= 0.19',
    'matplotlib >= 2.0',
    'seaborn >= 0.8'
    ],
)
