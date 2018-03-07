#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# class query external astronomical catalogs. Each instance of this class
# will be connected to (and search in) a specific catalog and collection.
# These catalogs have to be registered in a database collection of exposed
# (meaning they are ready to be used) catalogs.
#
# Once the instance is created it provides methods to query for:
#
#   - all the sources with a certain distance
#   - closest source at a given position
#   - binary search: return yes/no if anything is around the positon.
#
# The basic query functions are defined in the query_utils module.
#
# Author: M. Giomi (matteo.giomi@desy.de)

import pymongo
from astropy.table import Table
from healpy import nside2resol
from extcats.catquery_utils import searcharound_HEALPix, searcharound_9HEALPix, \
                            searcharound_2Dsphere, searcharound_RAW, get_closest


class CatalogQuery():
    """
        class to query external catalogs.
    """
    
    def __init__(self, cat_name, ra_key, dec_key, coll_name = "srcs", 
        dbclient = None, logger =  None):
        """
            Connect to the desired database and collection. Retrive information
            on how to query it and what's the catalog is about.
            
            Parameters:
            -----------
            
                cat_name: `str`
                    name of the catalog and database you want to query.
                
                ra[dec]_key: `str`
                    names of the coordinates in the catalog that will be used 
                    for the queries.
                
                coll_name: `str`
                    name of the catalog collection containing the sources.
                
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
        self.dbclient = dbclient
        if dbclient is None:
            self.dbclient = pymongo.MongoClient()
        self.logger.info("using mongo client at %s:%d"%(self.dbclient.address))
        
        # find database and collection
        if not cat_name in self.dbclient.database_names():
            raise KeyError("cannot find database %s in client."%(cat_name))
        self.cat_db = self.dbclient[cat_name]
        if not coll_name in self.cat_db.collection_names():
            raise KeyError("cannot find collection %s in database %s"%(coll_name, self.cat_db.name))
        self.src_coll = self.cat_db[coll_name]
        self.logger.info("connected to collection %s of database %s."%(coll_name, self.cat_db.name))
        self.logger.info("found %d documents in source collection %s."%(self.src_coll.count(), coll_name))
        
        # read metadata for the catalog
        if not "meta" in self.cat_db.collection_names():
            raise KeyError("cannot find metadata collection in database %s"%self.cat_db.name)
        
        # check for healpix and sphere2d support
        self.check_healpix()
        self.check_sphere2d()
        
        # check if query relevant keys are contained in the database
        self.ra_key = ra_key
        self.dec_key = dec_key
        important_keys  = [self.ra_key, self.dec_key]
        if self.has_hp:
            important_keys.append(self.hp_key)
        if self.has_2dsphere:
            important_keys.append(self.s2d_key)

        test_doc = self.src_coll.find_one()
        for k in important_keys:
            if not k in test_doc.keys():
                raise KeyError("key %s not found among document fields: %s"%(k, ", ".join(test_doc.keys())))
        
        # check which index does the collection actually have
        indexes = [k[:-2] for k in self.src_coll.index_information().keys()]
        self.logger.info("source collection has the following indexes: %s"%", ".join(indexes))
        if self.has_2dsphere and (not self.s2d_key in indexes):
            self.logger.warning("2dsphere key %s is not indexed."%(self.s2d_key))
            self.hp_index = False
        if self.has_hp and (not self.hp_key in indexes):
            self.logger.warning("2dsphere key %s is not indexed."%(self.hp_key))
            self.sphere2d_index = False


    def check_healpix(self):
        """
            reads metadata collection and figures out if the catalog support healpix.
            If there is just one healpix key, use it if indexed. If there are 
            more than one, assign the healpix key with the highest order that is indexed
        """
        
        hp_doc = None
        hp_docs = [doc for doc in self.cat_db["meta"].find({"type": "healpix"})]
        if len(hp_docs) == 0:
            self.logger.warning("no HEALPix description found in metadata for catalog %s"%self.cat_db.name)
        if len(hp_docs) == 1 and hp_docs[0]["is_indexed"]:
            hp_doc=hp_docs[0]
        else:
            self.logger.info("found %d HEALPix partitions for this catalog."%len(hp_docs))
            for i, _ in enumerate(hp_docs): self.logger.info("#%d %s"%(i,repr(_)))
            self.logger.info("Using the higher-order one of those that are indexed.")
            max_order = 0.
            for hpd in hp_docs:
                if hpd["is_indexed"] and hpd["order"] > max_order:
                    hp_doc = hpd
                    max_order = hpd["order"]

        if not hp_doc is None:
            self.has_hp = True
            self.hp_key, self.hp_order, self.hp_nest, self.hp_index = (
                hp_doc['key'], hp_doc['order'], hp_doc['nest'], hp_doc['is_indexed'] )
            self.hp_resol = nside2resol(2**self.hp_order, arcmin = True) * 60.
            self.logger.info("set HEALPIX partition of order %d with key %s. Nested: %s, Inexed: %s, Resolution [\"]: %.3f"%(
                self.hp_order, self.hp_key, str(self.hp_nest), str(self.hp_index), self.hp_resol))
        else:
            self.has_hp = False
            self.logger.warning("no indexed HEALPix partition found for this catalog.")


    def check_sphere2d(self):
        """
            reads metadata collection and figures out if the catalog support queries
            on spherical geometry. There should be at most one indexed 2dsphere-type
            of key (mongodb rule).
        """

        sphere2d_doc = [doc for doc in self.cat_db["meta"].find({"type": "sphere2d"})]
        if len(sphere2d_doc) == 0:
            self.has_2dsphere = False
            self.logger.warning("no 2d sphere key found in catalog %s"%self.cat_db.name)
        elif len(sphere2d_doc)==1 and sphere2d_doc[0]["is_indexed"]:
            # mongo collections can have at most one 2dsphere key indexed
            self.has_2dsphere = True
            self.sphere2d_index = True
            self.s2d_key = sphere2d_doc[0]["key"]
            self.logger.info("set 2dsphere key %s with format %s. Inexed: %s"%(
                self.s2d_key, sphere2d_doc[0]["pos_format"], self.sphere2d_index))
        else:
            self.logger.warning("mongo collections can have at most one 2dsphere key indexed.")


    def autoset_method(self):
        """
            decide the default method to use for queries based on the indexed keys
            found in the database. The idea is the following:
                
                - if there is an indexed geoJSON/'legacy pair' key use the 2dsphere methods.
                - elif there is an indexed healpix id use this one
                - elif there is an unindexed healpix use it (and raise Warning)
                - elif there is an unindexed geoJSON/'legacy pair' use it (and raise Warning)
                - else use the raw method (and raise Warning).
        """
        
        if self.sphere2d_index:
            self.default_method = "2dsphere"
        elif self.hp_index:
            self.default_method = "healpix"
        elif self.has_hp:
            self.logger.warning("no HEALPix / 2dsphere index found. Queries will be SUPER SLOW!")
            self.default_method = "healpix"
        elif self.has_2dsphere:
            self.logger.warning("no HEALPix / 2dsphere index found. Queries will be SUPER SLOW!")
            self.default_method = "2dsphere"
        else:
            self.logger.warning(
                "database doesn't even have any HEALPix / 2dsphere key. It will take ages to query.")
            self.default_method = "raw"
        self.logger.info("setting default search method to '%s'"%self.default_method)


    def findwithin_HEALPix(self, ra, dec, rs_arcsec, circular = True, find_one = False):
        """
            Returns sources in catalog contained in the the group of healpixels
            around the target coordinate that covers the search radius.
            
            NOTE that this query does not return stuff in a circle, but rather in 
            a square-like pattern of HEALpixels that covers the search radius.
            
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates of the target direction. They have to be in
                    the same format as those used to find the healpix id.
                
                rs_arcsec: `float`
                    maximum allowed distance to target position (in arcsec).
                
                circular: `bool`, default: True
                    if True, the results are skimmed removing sources outside of the
                    circular search radius.
                
                find_one: `bool`
                    if True the collection is searched with the find_one method returning
                    just the first result of the query. if False (default), the method
                    find is used, returning all matching documents.
                
            Returns:
            --------
                
                cptable: `astropy.table.Table`/None
                    astropy table of the catalog entry for the found counterpart. If
                    no counpterpart is found returns None.
        """
        
        if not self.has_hp:
            raise RuntimeError("catalog has no healpix index. Cannot use it to query.")
        
        return searcharound_HEALPix(
            ra = ra, dec = dec, rs_arcsec = rs_arcsec, src_coll = self.src_coll,
            hp_key = self.hp_key, hp_order = self.hp_order, 
            hp_nest = self.hp_nest, hp_resol = self.hp_resol, 
            circular = circular, ra_key = self.ra_key, dec_key = self.dec_key, find_one = find_one)


    def findwithin_9HEALPix(self, ra, dec, find_one = False):
        """
            Returns sources in catalog contained in the 9 healpixels
            around the target coordinate.
            
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates of the target direction. They have to be in
                    the same format as those used to find the healpix id.
                
                find_one: `bool`
                    if True the collection is searched with the find_one method returning
                    just the first result of the query. if False (default), the method
                    find is used, returning all matching documents.
                
            Returns:
            --------
                
                cptable: `astropy.table.Table`/None
                    astropy table of the catalog entry for the found counterpart. If
                    no counpterpart is found returns None.
        """
        
        if not self.has_hp:
            raise RuntimeError("catalog has no healpix index. Cannot use it to query.")
        self.logger.warning(
        "queries sources in a 9-pixel square of %.3f arcsec side around target"%(3*self.hp_resol))
        return searcharound_9HEALPix(
            ra = ra, dec = dec, src_coll = self.src_coll, hp_key = self.hp_key, 
            hp_order = self.hp_order, hp_nest = self.hp_nest, hp_resol = self.hp_resol, find_one = find_one)


    def findwithin_2Dsphere(self, ra, dec, rs_arcsec, find_one = False):
        """
            Returns sources in catalog within rs_arcsec from target position.
            
            This function uses the special mongodb queries on spherical geometry 
            using "$geoWithin" and "$centerSphere" operators. It requires the 
            catalog documents to have been assigned a geoJSON or 'legacy pair' 
            field of type 'Point' (see insert_example notebook).
            
            The Right Ascension (ra) should be within -180 and 180 degrees. 
            If not, it will be folded as ra = ra - 360 if ra > 180.
            
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates of the target direction. They have to be in
                    the same format as those used to build the 2dsphere object 
                    in the catalog documents.
                
                rs_arcsec: `float`
                    maximum allowed distance to target position (in arcsec).
                
                find_one: `bool`
                    if True the collection is searched with the find_one method returning
                    just the first result of the query. if False (default), the method
                    find is used, returning all matching documents.
                
            Returns:
            --------
                
                cptable: `astropy.table.Table`/None
                    astropy table of the catalog entry for the found counterpart. If
                    no counpterpart is found returns None.
        """
        if not self.has_2dsphere:
            raise RuntimeError("catalog has no geoJSON/legacy pair object. Cannot use it to query.")
        return searcharound_2Dsphere(
            ra = ra, dec = dec, rs_arcsec = rs_arcsec, 
            src_coll = self.src_coll, s2d_key = self.s2d_key, find_one = find_one)


    def findwithin_RAW(self, ra, dec, rs_arcsec, box_scale = 2., find_one = False):
        """
            Returns sources in catalog within rs_arcsec from target position.
            
            It first selects points within a box of radius 3 times larger than the
            search radius using $gte and $lte operators, then uses the $where expression
            to compute the angular distance of the sources in the box from the target.
            
            Use if you don't have any healpix in the catalog nor any 
            special coordinates that supports mongo 2dsphere queries.
            
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates (J2000, in degrees) of the target direction.
                
                rs_arcsec: `float`
                    maximum allowed distance to target position (in arcsec).
                
                box_scale: `float`
                    size of the square used to pre-search candidates, in units of rs_arcsec.
                
                find_one: `bool`
                    if True the collection is searched with the find_one method returning
                    just the first result of the query. if False (default), the method
                    find is used, returning all matching documents.
                
            Returns:
            --------
                
                cptable: `astropy.table.Table`/None
                    astropy table of the catalog entry for the found counterpart. If
                    no counpterpart is found returns None.
        """
        return searcharound_RAW(
            ra = ra, dec = dec, rs_arcsec = rs_arcsec, src_coll = self.src_coll,
            ra_key = self.ra_key, dec_key = self.dec_key, box_scale = box_scale, find_one = find_one)


    def findwithin(self, ra, dec, rs_arcsec, method = "healpix", **qfunc_args):
        """
            return all the catalog sources within rs_arcsec arcsecond from target
            positoin ra, dec (degrees).
        
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates (in degrees) of the target direction.
                
                rs_arcsec: `float`
                    maximum allowed distance to target position (in arcsec).
                
                method: `str`
                    how to query the database. If
                        - healpix: use the healpix index of the catalog (findwithin_HEALPix)
                        - 2dsphere: use mongodb searches in spherical geometry 
                        (findwithin_2Dsphere). 
                        - raw: run a javascript function to compute the angular distances
                        (use findwithin_RAW)
                
                qfunc_args: `dict`
                    additional arguments to be passed to the specific query function.
            
            Returns:
            --------
                
                cptable: `astropy.table.Table`/None
                    astropy table of the catalog entries for the sources in the search radius.
                    if no sources are found, returns None.
                
        """
        if method == "healpix":
            return self.findwithin_HEALPix(ra, dec, rs_arcsec, **qfunc_args)
        elif method == "2dsphere":
            return self.findwithin_2Dsphere(ra, dec, rs_arcsec, **qfunc_args)
        elif method == "raw":
            return self.findwithin_RAW(ra, dec, rs_arcsec, **qfunc_args)
        else:
            raise ValueError(
                "invalid method parameter: %s. Valid ones are 'healpix'/'2dsphere'/'raw'"%method)


    def findclosest(self, ra, dec, rs_arcsec, method = "healpix", **qfunc_args):
        """
            find the closest source to target ra and dec (degrees) and
            within rs_arcsec.
        
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates (in degrees) of the target direction.
                
                rs_arcsec: `float`
                    maximum allowed distance to target position (in arcsec).
                
                method: `str`
                    how to query the database. If
                        - healpix: use the healpix index of the catalog (findwithin_HEALPix)
                        - 2dsphere: use mongodb searches in spherical geometry 
                        (findwithin_2Dsphere). In this case a value for the sphere2d_key 
                        has to be passed to the function. 
                        - raw: run a javascript function to compute the angular distances
                        (use findwithin_RAW)
                
                qfunc_args: `dict`
                    additional arguments to be passed to the specific query function.
            
            Returns:
            --------
                
                This function returns a tuple (cptable, dist). If no macthes are found,
                the output is (None, None). Else:
                
                cptable: `astropy.table.Table`/None
                    astropy table of the catalog entry for the found counterpart.
                
                dist: `float`/None
                    angular distance (in arcsec) of the closest source to the target.
        """
        
        nearby = self.findwithin(
            ra = ra, dec = dec, rs_arcsec = rs_arcsec, method = method, **qfunc_args)
        if nearby is None:
            return None, None
        closest, min_dist = get_closest(
            ra =  ra, dec = dec, table = nearby, ra_key = self.ra_key, dec_key = self.dec_key)
        if min_dist>rs_arcsec:
            return None, None
        return closest, min_dist


    def binaryserach(self, ra, dec, rs_arcsec, method = "healpix", **qfunc_args):
        """
            Boolean query, returns True if any source has been found rs_arcsec
            arcseconds within target position.
        
            Parameters:
            -----------
                
                ra/dec: `float`
                    sky coordinates (in degrees) of the target direction.
                
                rs_arcsec: `float`
                    maximum allowed distance to target position (in arcsec).
                
                method: `str`
                    how to query the database. If
                        - healpix: use the healpix index of the catalog (findwithin_HEALPix)
                        - 2dsphere: use mongodb searches in spherical geometry 
                        (findwithin_2Dsphere). In this case a value for the sphere2d_key 
                        has to be passed to the function. 
                        - raw: run a javascript function to compute the angular distances
                        (use findwithin_RAW)
                
                qfunc_args: `dict`
                    additional arguments to be passed to the specific query function.
                
            Returns:
            --------
                
                isanything: `bool`
                    True if catalog sources are found within search radius.
        """
        
        nearby = self.findwithin(
            ra = ra, dec = dec, rs_arcsec = rs_arcsec, method = method, find_one = True, **qfunc_args)
        if nearby is None:
            return False
        else:
            return True


    def test_queries(self, query_type, method, rs_arcsec, npoints=1e4, points=None, rnd_seed = 42):
        """
            run test queries using a set of uniformly distributed points on 
            a sphere as targets. 
            
            Parameters:
            -----------
            
            query_type: `str`
                identifier for the type of query you want to test. Can either be
                'within' (for the findwithin query), 'closest' (for the findclosest), 
                or 'binary' (for the isanything query).
            
            method: `str`
                    how to query the database. If
                        - healpix: use the healpix index of the catalog (findwithin_HEALPix)
                        - 2dsphere: use mongodb searches in spherical geometry 
                        (findwithin_2Dsphere). In this case a value for the sphere2d_key 
                        has to be passed to the function. 
                        - raw: run a javascript function to compute the angular distances
                        (use findwithin_RAW)
            
            rs_arcsec: `float`
                    search radius in arcseconds.
            
            npoints: `int`
                number of MC generated points to be used for the macthing. 
            
            points: `None` or list/array-like
                use this argument to provide points to be used for the testing, 
                skipping the MC generation. points should be array like and have this
                structure: [[ra1, dec1], [ra2, dec2], ...]
            
            rnd_seed: `float` or None
                if not None, this seed will be passed to np.random.
            
            Returns:
            --------
            av_query_time: `float`
                average query time for npoints queries measured as (start-stop)/npoints
        """
        
        # method specif imports
        import time, tqdm, astropy
        
        # generate the random sample of points if none is given
        if points is None:
            from extcats.catquery_utils import random_point_sphere
            points = random_point_sphere(int(npoints), rnd_seed)
            self.logger.info("running test queries using %d random points"%(npoints))
        else:
            npoints = len(points)
            self.logger.info("running test queries using %d user defined points"%(npoints))
        
        # set up the query function
        if 'within' in query_type:
            qfunc = self.findwithin
        elif 'closest' in query_type:
            qfunc = self.findclosest
        elif 'binary' in query_type:
            qfunc = self.binaryserach
        else:
            raise ValueError(
            "illegal value for parameter 'query_type': %s. Allowed: 'within', 'closest', 'binary'."%
            query_type)
        self.logger.info("running %d test queries using function: %s and method: %s"%
            (npoints, qfunc.__name__, method))
        
        # measure query time
        tot_found = 0
        start = time.time()
        for pp in tqdm.tqdm(points):
            buff = qfunc(pp[0], pp[1], rs_arcsec, method)
            if ( type(buff) == tuple and buff == (None, None) ) or (buff is None) or (not buff):
                continue
            tot_found += 1
        end = time.time()
        total = end - start
        av_query_time = total / float(npoints)
        self.logger.info("Total document found in queries: %d"%tot_found)
        self.logger.info("Took %.2e sec for %d random queries. Average query time: %.3e sec\n"%
            (total, int(npoints), av_query_time))
    

    def userquery(self, qfilter, totable = True):
        """
            perform a user defined query on the database.
            
            Parameters:
            -----------
                
                qfilter: `str`
                    query expression
                
                totable: `bool`
                    weather to return the query results as astropy table or
                    raw documents (dictionaries)
            
            Returns:
            --------
                
                query results, either as astropy Table or list of raw documents.
                if no documents matches the query expression returns None.
        """
        
        qresults = [o for o in self.src_coll.find(qfilter)]
        if len(qresults)==0:
            return None
        else:
            if totable:
                return Table(qresults)
            else:
                return qresults



