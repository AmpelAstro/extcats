#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# download some test data to run example notebook
#
# Author: M. Giomi (matteo.giomi@desy.de)

import os

from urllib.request import urlretrieve
from shutil import unpack_archive



class MQC:
    
    def download(datadir = "./testdata/milliquas", 
        fname = "heasarc_milliquas.tdat.gz", overwrite = False):
        
        print ("downloading the million quasar catalog:")
        print ("ref: https://heasarc.gsfc.nasa.gov/W3Browse/all/milliquas.html")
        
        if not os.path.isdir(datadir):
            print ("creating testdata directory: %s"%datadir)
            os.makedirs(datadir)
            
        dest = os.path.join(datadir, fname)
        if (not os.path.isfile(dest)) and (not overwrite):
            urlretrieve("https://desycloud.desy.de/index.php/s/AsxkEJUJ8LrgJkd/download", dest)
            print ("succesfully downloaded the MQC to %s"%dest)
        else:
            print ("found test MQC data in %s. Use overwrite to re-download it."%dest)
        return dest

class PS1Small:
    
    def download(datadir = "./testdata/PS1DR1_test", overwrite = False):
        
        print ("downloading few PS1 DR1 MeanObject table files.")
        print ("ref: https://panstarrs.stsci.edu/")
        
        if not os.path.isdir(datadir):
            print ("creating testdata directory: %s"%datadir)
            os.makedirs(datadir)
        if not os.listdir(datadir) and (not overwrite):
            tmp = os.path.join(datadir, "tmp_PS1DR1_test.zip")
            urlretrieve("https://desycloud.desy.de/index.php/s/3AAJ30J7YZ7dH43/download", tmp)
            datadir_parent = os.path.abspath(os.path.join(datadir, os.pardir))
            unpack_archive(tmp, datadir_parent)
        else:
            print ("found test PS1DR1 data in %s. Use overwrite to re-download it."%datadir)
        return (datadir)
