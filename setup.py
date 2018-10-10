#!/usr/bin/env python
#
# This file is part of pyESVEP.
# Copyright 2018 Radoslaw Guzinski and contributors listed in the README.md file.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
from setuptools import setup

PROJECT_ROOT = os.path.dirname(__file__)

def read_file(filepath, root=PROJECT_ROOT):
    """
    Return the contents of the specified `filepath`.

    * `root` is the base path and it defaults to the `PROJECT_ROOT` directory.
    * `filepath` should be a relative path, starting from `root`.
    """
    try:
        # Python 2.x
        with open(os.path.join(root, filepath)) as fd:
            text = fd.read()
    except UnicodeDecodeError:
        # Python 3.x
        with open(os.path.join(root, filepath), encoding = "utf8") as fd:
            text = fd.read()    
    return text


LONG_DESCRIPTION = read_file("README.md")
SHORT_DESCRIPTION = "End-member-based Soil and Vegetation Energy Partitioning (ESVEP) model to estimate sensible and latent heat flux (evapotranspiration) from radiometric surface temperature data"
REQS = ['numpy>=1.10']

setup(
    name                  = "pyESVEP",
    packages              = ['pyESVEP'],
    dependency_links      = ['git+https://github.com/hectornieto/pyTSEB@master#egg=pyTSEB-0'],
    install_requires      = REQS,
    version               = "1.0",
    author                = "Radoslaw Guzinski",
    author_email          = "rmgu@dhigroup.com",
    maintainer            = "Radoslaw Guzinski",
    maintainer_email      = "rmgu@dhigroup.com",
    description           = SHORT_DESCRIPTION,
    license               = "GPL",
    url                   = "https://github.com/DHI-GRAS/py-esvep",
    long_description      = LONG_DESCRIPTION,
    classifiers           = [
        "Development Status :: Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Agricultural Science",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3"],
    keywords             = ['ESVEP','End-member-based Soil and Vegetation Energy Partitioning',
                            'Resistance Energy Balance', 'evapotranspiration', 'Remote Sensing',
                            'land surface energy fluxes'])
