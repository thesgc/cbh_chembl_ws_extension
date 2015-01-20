#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import cbh_chembl_ws_extension

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

version = cbh_chembl_ws_extension.__version__

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    print("You probably want to also tag the version now:")
    print("  git tag -a %s -m 'version %s'" % (version, version))
    print("  git push --tags")
    sys.exit()

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='cbh_chembl_ws_extension',
    version=version,
    description="""An extension to chembl web services to allow additonal functionality for write operations, permissions and collaboration""",
    long_description=readme + '\n\n' + history,
    author='Andrew Stretton',
    author_email='andrew.stretton@sgc.ox.ac.uk',
    url='https://github.com/strets123/cbh_chembl_ws_extension',
    packages=[
        'cbh_chembl_ws_extension',
    ],
    include_package_data=True,
    install_requires=[
        'django-grappelli==2.4.12',
        'xlsxwriter',
    ],
    license="MIT",
    zip_safe=False,
    keywords='cbh_chembl_ws_extension',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
    ],
)
