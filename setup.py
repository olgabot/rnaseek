#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    'pandas',
]

test_requirements = [
    # TODO: put package test requirements here
]

scripts = ['rnaseek/combine_sailfish_output.py',
           'rnaseek/combine_miso_output.py']

setup(
    name='rnaseek',
    version='0.1.0',
    description='Library for parsing, aggregating, and overall dealing with '
                'outputs from RNA-sequencing data files',
    long_description=readme + '\n\n' + history,
    author='Olga Botvinnik',
    author_email='olga.botvinnik@gmail.com',
    url='https://github.com/olgabot/rnaseek',
    packages=[
        'rnaseek',
    ],
    scripts=scripts,
    package_dir={'rnaseek':
                 'rnaseek'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='rnaseek',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
