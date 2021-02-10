#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# setup script for the extcats package.
#
# Author: M. Giomi (matteo.giomi@desy.de)

from setuptools import setup

setup(
    name="extcats",
    version="2.4",
    description="Tools to organize and query astronomical catalogs",
    author="Matteo Giomi",
    author_email="matteo.giomi@desy.de",
    packages=["extcats"],
    url="https://github.com/AmpelProject/extcats",
    download_url="https://github.com/AmpelProject/extcats/archive/2.4.tar.gz",
    install_requires=["pymongo>=3.7", "healpy", "astropy", "pandas", "tqdm"],
    extras_require={"testing": ["pytest", "pytest-cov"]},
)
