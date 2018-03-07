#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# setup script for the extcats package.
#
# Author: M. Giomi (matteo.giomi@desy.de)

import os
from urllib.request import urlretrieve
from shutil import unpack_archive
from setuptools import setup

setup(
    name='extcats',
    version='1.4',
    description='Tools to organize and query astronomical catalogs',
    author='Matteo Giomi',
    author_email='matteo.giomi@desy.de',
    packages=['extcats'],
    url = 'https://github.com/MatteoGiomi/extcats',
    download_url = 'https://github.com/MatteoGiomi/extcats/archive/1.3.tar.gz',
    install_requires=['pymongo', 'healpy', 'astropy', 'pandas', 'tqdm'],
    )

# ------- download test data ------- #
# the million quasar catalog
mqc_testdata_dir = "./testdata/milliquas"
if not os.path.isdir(mqc_testdata_dir):
    print ("creating testdata directory: %s"%mqc_testdata_dir)
    os.makedirs(mqc_testdata_dir)
dest = os.path.join(mqc_testdata_dir, "heasarc_milliquas.tdat.gz")
if not os.path.isfile(dest):
    print ("downloading the million quasar catalog \
(MQC, https://heasarc.gsfc.nasa.gov/W3Browse/all/milliquas.html)")
    urlretrieve("https://desycloud.desy.de/index.php/s/AsxkEJUJ8LrgJkd/download", dest)

# a few small files from Panstarrs
ps1_testdata_dir = "./testdata/PS1DR1_test"
if not os.path.isdir(ps1_testdata_dir):
    print ("creating testdata directory: %s"%ps1_testdata_dir)
    os.makedirs(ps1_testdata_dir)
if not os.listdir(ps1_testdata_dir):
    print ("downloading few PS1 DR1 MeanObject table files (https://panstarrs.stsci.edu/)")
    tmp = os.path.join(ps1_testdata_dir, "tmp_PS1DR1_test.zip")
    urlretrieve("https://desycloud.desy.de/index.php/s/3AAJ30J7YZ7dH43/download", tmp)
    unpack_archive(tmp, "./testdata")
