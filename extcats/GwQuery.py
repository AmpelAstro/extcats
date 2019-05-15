#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# class to match coordinates to GW contour shapes.
#
# Author: M. Giomi (matteo.giomi@desy.de)

import pymongo

class GwQuery():
    """
         class to match coordinates to GW contour shapes.
    """
    
    def __init__(self, db_name, coll_name = "srcs", 
        dbclient = None, logger =  None):
        """
            Connect to the desired database and collection. Retrive information
            on how to query it and what's inside.
            
            Parameters:
            -----------
            
                db_name: `str`
                    name of database where the GW contours are saved.
                
                coll_name: `str`
                    name of the catalog collection containing the contours.
                
                dbclient: `pymongo.mongo_client.MongoClient`
                    mongodb client where. If None, it will use the default one.
                
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
        
        # connect to mongo client
        self.dbclient = pymongo.MongoClient() if dbclient is None else dbclient
        self.logger.debug("using mongo client at %s:%d"%(self.dbclient.address))
        
        # find database and collection
        if not db_name in self.dbclient.database_names():
            raise KeyError("cannot find database %s in client."%(cat_name))
        self.gw_db = self.dbclient[db_name]
        if not coll_name in self.gw_db.collection_names():
            raise KeyError("cannot find collection %s in database %s"%(coll_name, self.gw_db.name))
        self.gw_coll = self.gw_db[coll_name]
        self.logger.info("connected to collection %s of database %s. %s documents in collection."%
            (coll_name, self.gw_db.name, self.gw_coll.count()))
        
        # read metadata for the catalog
        if not "meta" in self.gw_db.collection_names():
            raise KeyError("cannot find metadata collection in database %s"%self.gw_db.name)
        
        # check the metadata for polygons keys
        self.check_polygon()


    def set_polygon_key(self, poly_key):
        """
            set polygon key.
        """
        self.poly_key = poly_key


    def check_polygon(self):
        """
            read the metadata collection to find the name of the key
            with the geoJSON elements representing the contours.
        """
        geo_doc = [doc for doc in self.gw_db["meta"].find({"type": "geoJSON.Polygon"})]
        if len(geo_doc) == 1:
            self.poly_key = geo_doc[0]['key']
            self.logger.debug("found geoJSON.Polygon key %s in metadata collection"%(self.poly_key))
        elif len(geo_doc) == 0:
            self.poly_key = None
            self.logger.warning("no geoJSON.Polygon key found in metadata collection")
        else:
            self.logger.warning("more than one geoJSON.Polygon key in metadata collection. Use set_polygon_key().")
            self.poly_key = None


    def query_contours(self, ra, dec, radius=None, **kwargs):
        """
            query the database for all contours containing the specified point.
            Optinally includes radius to account for positional uncertainties on
            ra, dec position.
            
            Parameters:
            -----------
            
                ra/dec: `float`
                    sky coordinates (in degrees) of the target direction.
                
                radius: `float`
                    uncertainty on ra, dec position, in arcseconds. 
                
                **kwargs:
                    to be passed to pymongo.Collection.find.
            
            Returns:
            --------
                
                list of documents matching your query.
        """
        
        if self.poly_key is None:
            raise ValueError("geoJSON.Polygon key not set. Use set_polygon_key().")
        
        # fold coordinates: geoJSON longitude goes from -180 to 180
        tdec = dec
        tra = ra if ra < 180 else (ra - 360)
        
        # buid the positional query
        query = {
                self.poly_key: {
                    '$geoIntersects': {
                        '$geometry': {'type': 'Point', 'coordinates': [tra, tdec]}
                    }
                }
            }
        
        # issue query and return the results
        return list(self.gw_coll.find(query, **kwargs))



