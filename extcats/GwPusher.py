#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# class to insert probability contours from GW maps into mongod
#
# Author: M. Giomi (matteo.giomi@desy.de)


import os, json, subprocess
import geojson
import copy
import numpy as np
from astropy.io import fits
from pymongo import MongoClient, GEOSPHERE

def map_to_feature_collection(gw_localization_map, confidence_levels=0.1):
    """
        create geoJSON FeatureCollection with the contours of
        GW position maps at given confidence level(s).
        
        Parameters:
        -----------
            
            gw_localization_map: `str`
                path to fits file containing the pdf of GW event location.
            
            confidence_levels: `float` or `array-like`
                value(s) of the confidence levels for which to compute the contours
        
        Returns:
        --------
            
            geojson.FeatureCollection object with the contours.
    """
    
    # parse confidence level(s)
    try:
        cls = " ".join([str(cl) for cl in confidence_levels])
    except TypeError:
        cls = str(confidence_levels)
    
    # execute command and obtain output
    cmd = subprocess.Popen(
                'ligo-skymap-contour %s --contour %s'%(gw_localization_map, cls),
                shell=True, 
                stdout=subprocess.PIPE
            )
    json_str, _ = cmd.communicate()
    
    # format into geojson, check and return
    fc = geojson.FeatureCollection(geojson.loads(json_str)['features'])
    if not fc.is_valid:
        raise ValueError('input string %s is not valid geoJSON. Errors: %s'%
                        (json_str, repr(fc.errors())))
    return fc


def fold_coords(feature, inplace=False):
    """
        fold longitude (R.A.) from 0-360 to -180 to 180.
        
        Parameters:
        -----------
        
            feature: `geojson.Feature`
                object for which the coordinates will be folded.
            
            inplace: `bool`
                if True, the passed object will be modified. Else, we work 
                on a copy.
        
        Returns:
        --------
            
            new geojson.Feature with folde coordinates if inplace is False, else None.
    """
    
    def fold(cp):
        if cp[0] > 180: cp[0] = cp[0] - 360
        return cp
    
    if inplace:
        geojson.utils.map_tuples(fold, feature)
        return None
    else:
        feature_copy = copy.deepcopy(feature)
        geojson.utils.map_tuples(fold, feature_copy)
        return feature_copy

def to_poly(coords):
    """
        cast list of coordinates to geojson.Polygon.
    """
    pp = geojson.Polygon(coords)
    if not pp.is_valid:
        raise ValueError("cannot convert coordinates %s to valid geojson.Polygon. Errors: %s"%
                         (repr(my_coords), pp.errors()))
    return pp


def feature_to_polygons(feature):
    """
        given a geoJSON Feature object, convert it into Polygon(s). 
        A Feature object represents a spatially bounded thing. 
        In the GW case, these Features contains MultiLineStrings. 
        
        We build the polygons from the coordinates of each of those.
        
        Parameters:
        -----------
        
            feature: `geojson.Feature`
                object for which the coordinates will be folded.
        
        Returns:
        --------
        
            list of geojson.Polygon objects
    """
    
    # first fold the coordinates, all in one go.
    # perform the operation on a copy of the feature.
    folded = fold_coords(feature)
    
    # now convert to multiline string
    mls = geojson.MultiLineString(folded['geometry'], validate=True)
    if not mls.is_valid:
        raise ValueError("cannot convert feature geometry %s to valid geojson.MultiLineString"%
                         (repr(feature.coordinates)))
     
    # convert each coordinate list of the multilinestring to polygons and return
    return [to_poly([cc]) for cc in mls['coordinates']]


class GwPusher():
    """
        Class to insert probability countours from GW localization maps into mongod
    """

    def __init__(self, dbname, coll_name, dbclient=None, logger=None):
        """
            Connect to the database where you want to ingest the GW contours.

            Parameters:
            -----------
                
                dbclient: `pymongo.mongo_client.MongoClient`
                    mongodb client. If None, it will use the default one.
                
                dbname: `str` or None
                    name of the database to used.
                
                coll_name: `str`
                    name of the collection. Default is 'gw_contours'.
                
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
        
        # connect to databse and collection
        self.dbclient = MongoClient() if dbclient is None else dbclient
        self.dbname = dbname
        self.coll_name = coll_name
        self.db = self.dbclient[dbname]
        self.coll = self.db[coll_name]
        
        # be nice and informative
        self.logger.info("connected to collection %s of database %s"%
            (self.dbname, self.coll_name))
        self.logger.debug("using mongo client at %s:%d"%(self.dbclient.address))
        self.logger.debug("database stats: %s"%json.dumps(self.db.command("dbstats"), indent=2))


    def insert_gw_contours(self, gw_localization_map, confidence_levels, gw_id=None, 
        geom_key="poly", overwrite_coll=False, **kwargs):
        """
            extract contours from GW localization map for the specified confidence
            levels and create JSON documents with polygon geometry to be inserted 
            in the database.

            Parameters:
            -----------
                
                gw_localization_map: `str`
                    path to fits file containing the pdf of GW event location.
            
                confidence_levels: `float` or `array-like`
                    value(s) of the confidence levels for which to compute the contours
                
                gw_id: `str`
                    ID of the GW event. If None, it will be read from the fits file.
                
                geom_key: `str`
                    name of the document key conatinig the geoJSON polygon. This is the 
                    key that will be used in the search.
                
                overwrite_coll: `bool`
                    if True and another collection with the given name exists, it will
                    be dropped and replaced.
                
                **kwargs:
                    to be passed to pymongo.Collection.insert_many
                
            Returns:
            --------
                docs: `list` of JSON documents extracted from the fits file.
        """
        
        self.logger.info("inserting contours at %s CL for %s in collection %s."%
            (gw_localization_map, repr(confidence_levels), ".".join([self.db.name, self.coll.name])))

        # get the ID of the map from the header
        if gw_id is None:
            gw_id = str(fits.getheader(gw_localization_map, 1)['OBJECT'])
        self.logger.info("using GW event ID: %s"%repr(gw_id))
        
        # compute the contours at the desired confidence levels
        fc = map_to_feature_collection(gw_localization_map, confidence_levels)
        
        # set the geom key so that we know what to use for the query / index
        self.geom_key = geom_key
        
        # create a document for all the features in the feature collection
        docs = []
        for feature in fc.features:
            
            # extract confidence level of the feature
            cl = feature['properties']['credible_level']
            
            # extract geometry as a Polygon(s)
            polygons = feature_to_polygons(feature)
            
            for poly in polygons:
                # create the docs and add it to the list
                buff = {
                    'gwID': gw_id,
                    'file': gw_localization_map,
                    'CL': cl,
                    self.geom_key: {
                                    "type": "Polygon",
                                    "coordinates": poly.coordinates
                                }
                    }
                docs.append(buff)
        self.logger.debug("extracted %d documents from GW map %s"%(len(docs), gw_localization_map))
        
        # now the insertion part: drop if you care
        if self.coll_name in self.db.collection_names() and overwrite_coll:
            self.logger.warning("overwrite_coll asserted: collection %s will be dropped."%self.coll_name)
            self.db.drop_collection(self.coll_name)
        
        # build the indexes
        self.coll.create_index([(self.geom_key, GEOSPHERE)])
        #TODO: add index for GW ID field.
        
        # insert the docs
        res = self.coll.insert_many(docs, **kwargs)
        self.logger.info("inserted %d documents (out of %d) into the database."%
            (len(res.inserted_ids), len(docs)))
        
        # add metadata fpr polygon key
        self.polygon_meta()
        
        # returned inserted documents
        return docs


    def polygon_meta(self):
        """
            Add metadata identifying the key with the Polygon coordinates.
            
                {
                    'key' : self.geom_key
                    'type': 'geoJSON.Polygon'
                }
        """
        
        self.logger.info("adding metadata for geoJSON Polygon field named: %s"%self.geom_key)
        self.db['meta'].update(
            {'_id': self.geom_key},
            {
              '_id' : self.geom_key,
              'key' : self.geom_key,
              'type': 'geoJSON.Polygon'
            },
             upsert = True)


    def science_meta(self, contact, email, description, reference):
        """
            Add metadata associated to the catalog as a tool for science.

            Paramters:
            ----------

                contact/email: `str`
                    name and email of contact person.

                description: `str`
                    brief description of the catalog content

                refernce: `str`
                    publication, web address documenting the catalog.
        """

        self.logger.info("adding science related metadata to database.")
        self.db['meta'].update(
            {'_id': 'science'},
            {
              '_id' : 'science',
              'contact': contact,
              'email': email,
              'description': description,
              'ref': reference,
            },
            upsert = True)

