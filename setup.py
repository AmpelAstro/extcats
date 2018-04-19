#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# setup script for the extcats package.
#
# Author: M. Giomi (matteo.giomi@desy.de)

from setuptools import setup

setup(
    name='extcats',
    version='1.8',
    description='Tools to organize and query astronomical catalogs',
    author='Matteo Giomi',
    author_email='matteo.giomi@desy.de',
    packages=['extcats'],
    url = 'https://github.com/MatteoGiomi/extcats',
    download_url = 'https://github.com/MatteoGiomi/extcats/archive/1.3.tar.gz',
    install_requires=['pymongo', 'healpy', 'astropy', 'pandas', 'tqdm'],
    )
