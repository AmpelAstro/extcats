#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# class to manage the interfacte between esternal catatalog databases
#
# Author: M. Giomi (matteo.giomi@desy.de)

import pymongo
import json

class CatalogManager():
    
    def __init__(self, exposed_cats_dbclient =  None, 
        exposed_cats_db = "exposed_catalogs", exposed_cats_coll = "cats", logger = None):
        """
            Paramaters:
            -----------
                
                dbclient: `pymongo.mongo_client.MongoClient` or None
                    mongodb client. If None, it will use the default one.
                
                exposed_cats_db: `str` or pymongo.database.Database
                    database where the exposed catalog collection is kept. If str
                    this database is supposed to be at the same client as the one
                    used for pushing. 
                
                exposed_cats_coll: `str`
                    name of the collection of the exposed_cats_db containing the
                    database names. 
                
                logger: `self.logger.Logger`:
                    logger for the class. If None, a default one will be created.
        """
        
        # init the logger
        if logger is None:
            import logging
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger
        
        # connect to database
        if exposed_cats_dbclient is None:
            exposed_cats_dbclient = pymongo.MongoClient()
        self.dbclient = exposed_cats_dbclient
        self.expcats_db = self.dbclient[exposed_cats_db]
        self.logger.info("connected to database %s"%exposed_cats_db)        
        self.logger.info("using mongo client at %s:%d"%self.dbclient.address)
        
        # connect to collection
        self.expcats_coll = self.expcats_db[exposed_cats_coll]
        self.logger.info("connected to collection %s. %d documents inside."%(
            exposed_cats_coll, self.expcats_coll.count()))


    def listcatalogs(self):
        """
            list all the catalogs in the exposed catalog collection.
        """
        
        self.logger.info("exposed cats collection contains: \n%s"%
            json.dumps([o for o in self.expcats_coll.find()], indent=2))


    def expose(self, cat_name):
        """
            add a given catalog to the exposed cats collection. Only catalogs 
            that are expose will be accessible to the query routines. 
            
            Paramaters:
            -----------
            
                catname: `str`
                    name of the catalog to expose.
        """
        
        if not self.isexposed(cat_name):
            self.expcats_coll.insert_one({"_id": cat_name })
            self.logger.info("added catalog %s to exposed_catalogs database:"%cat_name)
        else:
            self.logger.warning("catalog %s was already exposed."%cat_name)
        self.listcatalogs()


    def hide(self, cat_name):
        """
            remove the given catalog from the exposed cats collection.
            
            Paramaters:
            -----------
            
                catname: `str`
                    name of the catalog to hide.
        """
        
        if self.isexposed(cat_name):
            self.expcats_coll.delete_one({"_id": cat_name })
            self.logger.info("removed catalog %s to exposed_catalogs database:"%cat_name)
        else:
            self.logger.warning("catalog %s was not exposed."%cat_name)
        self.listcatalogs()


    def isexposed(self, cat_name):
        """
            find out if a catalog is present in the exposed cat collection. Returns
            a boolean.
            
            Paramaters:
            -----------
            
                catname: `str`
                    name of the catalog to expose.
        """
        
        return bool(self.expcats_coll.find_one({"_id": cat_name }))
