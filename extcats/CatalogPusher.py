#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# class to insert static esternal catalogs into a mongodb.
#
# Author: M. Giomi (matteo.giomi@desy.de)

import os, glob, time, json, pymongo
import pandas as pd
import numpy as np
from inspect import getsourcelines


class CatalogPusher():
    """
        Class to insert static external catalogs into a mongodb.
    """


    def __init__(self, catalog_name, data_source, file_type = None, recursive = True, logger = None):
        """
            Parameters:
            -----------
            
            catalog_name: `str`
                short name of the catalog.
            
            data_source: `str` or `list`
                specifies where the raw files fo the catalog are located. Either
                single string or list of strings are accepted. If data_source is 
                a string, it can either be the path of a single file, or a directory.
                In case it is a directory, all the files inside that directory (
                and its subdirectories) are considered. If data_source is a list 
                of directories, each directory will be searched for files.
            
            file_type: `str` or None
                extension (e.g. .csv, .csv.gz, .fits, ecc) of the raw files to
                be ingested. 
            
            recursive: `bool`
                if True and data_source contains directories, then they will be 
                searched recursively for files.
            
            logger: `self.logger.Logger`:
                logger for the class. If None, a default one will be created.
        """
        
        # and man gave names
        self.cat_name = catalog_name
        
        # init the logger
        if logger is None:
            import logging
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger
        
        # find the raw files for this catalog. If a single string has been
        # passed, cast it into list so we treat the two cases the same way
        if type(data_source) == str:
            data_source = [data_source]
        elif type(data_source) == list:
            pass
        else:
            raise OSError("data_source param must be string or lits, got %s instead"%type(data_source))
        
        # check files / search for directories
        self.raw_files = []
        for path in data_source:
            if os.path.isfile(path):
                self.raw_files.append(path)
            elif os.path.isdir(path):
                file_filter = ("/**/*"*recursive + "/*"*(not recursive))
                if not file_type is None:
                    file_filter += file_type
                files = glob.glob(path+file_filter, recursive = recursive)
                self.raw_files.extend(files)
            else:
                raise OSError("invalid path to raw files: %s"%path)
        if len(self.raw_files) == 0:
            self.logger.warning("no data found for catalog %s in source: %s"%(catalog_name, data_source))
        else:
            self.logger.info("found %d files for catalog %s in data source: %s"%(
                    len(self.raw_files), self.cat_name, data_source))
        
        # check them from existence and consistency
        self.raw_files.sort()
        self.check_raw_files()


    def check_raw_files(self):
        """
            perform checks on the raw files for the catalog. They have to exist, 
            and to have all the same extension. 
        """
        
        self.logger.info("checking raw files for existence and consistency..")
        
        for rawfile in self.raw_files:
            if not os.path.isfile(rawfile):
                raise OSError("invalid file path: %s"%rawfile)
        
        types_set = set([os.path.splitext(rawf)[-1] for rawf in self.raw_files])
        if len(types_set) > 1:
            raise TypeError("inconsistent file types %s in raw_files: \n%s.\nUse file_type parameter in constructor."%
                    (", ".join(types_set), "\n".join(self.raw_files)))
        self.logger.info("all files exists and have consistent type.")
    
    def nfiles(self):
        """
            return the number of raw_files in the object.
        """
        return len(self.raw_files)


    def file_groups(self, group_size = 4):
        """
            returns indexes that can be used to divide the list
            of raw_files into groups of (at max) group_size files. 
            The last group can be smaller if len(self.raw_files)%group_size != 0.
        
            Parameters:
            -----------
            
                group_size: `int`
                    size of file groups, equal to x[1]-x[0] where x is an element of the
                    output list.
            
            Returns:
            --------
                    
                    output Nx2 list, each element of the list represent a valid slice for the
                    self.raw_files file list.
        """
        
        return [[x, x+group_size] for x in range(0, len(self.raw_files), group_size)]


    def assign_file_reader(self, reader_func, read_chunks, **reader_args):
        """
            define the function used to read the single files to be ingested in
            the catalog. 
            
            Parameters:
            -----------
            
                readerfunc: `callable`
                    function that will be used to read the raw files and insert them
                    into the database. Supported functions are pandas.DataFrame.read_* 
                    or, for fits files, astropy.table.Table
                
                read_chunks: `bool`
                    weather or not the supplied reader support reading raw files in chunks, 
                    for example if the reader is one of the pandas read_table[csv/fwf/ecc]
                    functions and the chunksize option is given.
                    
                reader_args: 
                    valid arguments for the reader_func.
        """
        
        # check if given keyword arguments are accepted by the function
        for rdr_arg in reader_args.keys():
            if not rdr_arg in reader_func.__code__.co_varnames:
                self.logger.warning(
                "argument %s not found among those accepted by %s"%(rdr_arg, reader_func.__name__))
        self.file_reader = reader_func
        self.file_reader_inchunks = read_chunks
        self.file_reader_args = reader_args
        self.logger.info("file reader %s assigned to pusher."%(self.file_reader.__name__))


    def assign_dict_modifier(self, modifier_func):
            """
                define the function used to format the documents (dictionaries) before
                inserting them into the database.
                
                Parameters:
                -----------
                
                    modifier_func: `callable`
                        function that act on the single dictionaries and format them for
                        db ingestion. This function must accept a dictionary and return a 
                        dictionary.
            """
            self.dict_modifier = modifier_func
            self.logger.info("source document modifer %s assigned to pusher."%
                (self.dict_modifier.__name__))


    def insert_file_todb(self, raw_file, db_coll, dry = False):
        """
            parse a raw file into a list of dictionaries that can be inserted in
            the database, and apply the dict modifier to each element. 
            
            Parameters:
            -----------
                raw_file: `str`
                    path to a file with catalog entries.
                
                db_coll: `pymongo.collection.Collection`
                    databse collection to be fed.
                
                dry: `bool`
                    if True, raw files are read and parsed into documents, but 
                    the database won't be filled.
            
            Returns:
            --------
                    list of the names of the processed file.
        """
        
        # check if the reader has been defined # TODO: try to guess from file extension and go.
        if not hasattr(self, "file_reader"):
            raise AttributeError(
            "no file_reader has been defined for this object. Use assign_file_reader first.")
        if not hasattr(self, "dict_modifier"):
            raise AttributeError(
            "no dict_modifier has been defined for this object. Use assign_dict_modifier first.")
        
        self.logger.info("inserting %s in collection %s."%
            (raw_file, ".".join([db_coll.database.name, db_coll.name])))
        
        
        # read raw file, convert into documents, and insert to DB.
        def convert_and_push(data, dry = dry):
            """
                helper function to treat the two cases (chunked/unchunked) the 
                same way. Convert data into dictionaries (use pandas default 
                functionality) if it's an astropy Table convert it into 
                pandas first. Returns the number of document inserted.
                
            """
            if hasattr(data, "to_dict"):
                raw_docs = data.to_dict("records")
            elif hasattr(data, "to_pandas"):
                raw_docs = data.to_pandas().to_dict("records")
            else:
                raise NotImplementedError(
                "only astropy Table or pandas DataFrame/TextReader are supported return types of file_reader.")
            docs = [self.dict_modifier(dd) for dd in raw_docs]
            if not dry:
                db_coll.insert_many(docs)
            return len(docs)
        
        # now read, parse, and fill.
        tot_docs, start = 0, time.time()
        if not self.file_reader_inchunks:
            data = self.file_reader(raw_file, **self.file_reader_args)
            tot_docs += convert_and_push(data)
        else:
            for data in self.file_reader(raw_file, **self.file_reader_args):
                tot_docs += convert_and_push(data)
        end = time.time()
        if not dry:
            self.logger.info("inserted %d documents in %.2e seconds"%(tot_docs, end-start))
        else:
            self.logger.info("DRY RUN: parsed %d documents in %.2e seconds"%(tot_docs, end-start))
        return raw_file


    def push_to_db(self, dbclient = None, dbname = None, coll_name = "srcs", 
        index_on = None, index_args = None, overwrite_coll = False,
        append_to_coll = True, dry = False, filerange = None):
        """
            insert all the files into the specified database.
            
            Parameters:
            -----------
            
                dbclient: `pymongo.mongo_client.MongoClient`
                    mongodb client. If None, it will use the default one.
                
                dbname: `str` or None
                    name of the database, if None the name of the catalog (self.cat_name)
                    will be used.
                
                coll_name: `str`
                    name of the collection. Default is 'srcs'.
                
                index_on: `str`, `list` or None
                    document key (keys if list is provided) for which an index will be built. 
                    If None, no index will be created for the collection.
                
                index_args: None or `list`
                    additional arguments to be passed to pymongo collection create_index
                    method. The list should have the same length as the number of indexes
                    to be created (index_args).
                
                overwrite_coll: `bool`
                    if True and another collection with the given name exists, it will
                    be dropped and replaced.
                
                append_to_coll: `bool`
                    if True and another collection with the given name exists, new documents
                    will be added to the collection. Else, nothing will be done.
                
                dry: `bool`
                    parameter is passed to insert_file_todb. If True, raw files are 
                    read and parsed into documents, but the database won't be filled.
                
        """
        
        # connect to databse and collection
        if dbclient is None:
            dbclient = pymongo.MongoClient()
        self.dbclient = dbclient
        self.logger.info("using mongo client at %s:%d"%(self.dbclient.address))
        
        if dbname is None:
            dbname = self.cat_name
        db = dbclient[dbname]
        self.dbname = dbname
        self.db = db
        self.logger.info("connecting to database %s. Here some stats:"%(db.name))
        self.logger.info(json.dumps(db.command("dbstats"), indent=2))
        
        if coll_name in db.collection_names() and overwrite_coll:
            self.logger.warning("overwrite_coll asserted: collection %s will be dropped."%coll_name)
            db.drop_collection(coll_name)
        coll = db[coll_name]
        self.coll = coll
        self.coll_name = coll_name
        if coll_name in db.collection_names() and (not append_to_coll):
            self.logger.warning("collection %s already exists and append_to_coll is False. Nothing to do."%coll_name)
            return 
        
        # create the indexes
        if not index_on is None:
            if type(index_on) == str:
                index_list = [index_on]
            elif type(index_on) == list:
                index_list = index_on
            for i_ind, index in enumerate(index_list):
                if index_args is None:
                    coll.create_index(index)
                else:
                    coll.create_index(index, **index_args[i_ind])
        self.logger.info("collection has the following indexes: %s"%
            ", ".join(coll.index_information().keys()))
        
        # push the files into the database, either all of them, or some with
        # specified index range.
        start = time.time()
        files = self.raw_files
        if (not filerange is None) and len(self.raw_files)>1:
            files = self.raw_files[filerange[0]:filerange[1]]
        for rawf in files:
            self.insert_file_todb(rawf, coll, dry)
        end = time.time()
        self.logger.info("done inserting catalog %s in collection %s.%s. Took %.2e seconds"%
            (self.cat_name, dbname, coll_name, end-start))


    def info(self):
        """
            print out information on the ingetsed database, such as the size of 
            the collection, its indexes and sizes, the format of the inserted 
            documents, and so on. This inforamtion is stored in the self.infodict
            attribute of the object.
        """
        
        self.logger.info("here are some statis on database %s and collection %s"%(self.dbname, self.coll_name))
        self.infodict =  {
            "cat_name": self.cat_name,
            "io_info":  {
                "raw_files": self.raw_files,
                "reader":   {
                    "name": self.file_reader.__name__,
                    "args": self.file_reader_args
                            },
                "modifier": {
                    "name": self.dict_modifier.__name__,
                    "code": getsourcelines(self.dict_modifier)
                            }
                        },
            "dbinfo":   {
                "client_address": self.dbclient.address,
                "dbstats": self.db.command("dbstats"),
                "src_coll_stats": self.db.command("collstats", self.coll_name)
                        }
                    }
        print(json.dumps(self.infodict, indent=2))


    def run_test(self, query_func, npoints=1e6, points=None):
        """
            test the database cross matching it with a set of fake points on the
            sphere.
            
            Parameters:
            -----------
            
                query_func: callable
                    function that is used to query the catalog. It must accept ra and
                    dec, and the collecion as the arguments:
                        objs = query_func(ra, dec)
                
                npoints: `int`
                    number of MC generated points to be used for the macthing. 
                
                points: `None` or list/array-like
                    use this argument to provide points to be used for the testing, 
                    skipping the MC generation. points should be array like and have this
                    structure: [[ra1, dec1], [ra2, dec2], ...]
            
            Returns:
            --------
                av_query_time: `float`
                    average query time for npoints queries.
        """
        
        # method specific import
        import tqdm
        
        # generate the random sample of points
        if points is None:
            from extcats.catquery_utils import random_point_sphere
            points = random_point_sphere(int(npoints))
            self.logger.info("running test queries using %d random points"%(npoints))
        else:
            npoints = len(points)
            self.logger.info("running test queries using %d user defined points"%(npoints))
        
        # measure query time
        tot_found = 0
        start = time.time()
        for pp in tqdm.tqdm(points):
            buff = query_func(ra = pp[0], dec = pp[1], coll = self.coll)
            if not buff is None:
                tot_found += 1
        end = time.time()
        total = end - start
        av_query_time = total / float(npoints)
        self.logger.info("Total document found for query: %d"%tot_found)
        self.logger.info("Took %.2e sec for %d random queries. Average query time: %.3e sec"%
            (total, int(npoints), av_query_time))


    def healpix_meta(self, healpix_id_key, order, is_indexed, nest = True):
        """
            Add metadata describing the healpix grid used to partition the sources.
            This info will be arranged in the metdata collection as a separate document
            with the structure:
            
                { 
                    'key' : healpix_id_key,
                    'order' : order
                    'nest': nest
                    'is_indexed': is_indexed
                }
            
            If multiple healpix paritions are available for the catalog, the healpix_id_key
            will be used to tag the documents describing each of them.
            
            Paramters:
            ----------
                
                helpix_id_key: `str`
                    name of healpix key in the documents.
                
                order: `int`
                    order of the healpix grid, so that nside = 2**order.
                
                is_indexed: `bool`
                    weather or not the srcs collection has been indexed using this key.
                
                nest: `bool` 
                    weather or not healpix grid has nested geometry. Default is True.
        """
        
        self.logger.info("adding metadata for HEALPix key named: %s"%healpix_id_key)
        self.db['meta'].update(
            {'_id': healpix_id_key},
            { 
              '_id' : healpix_id_key,
              'key' : healpix_id_key,
              'order' : order,
              'nest': nest,
              'is_indexed': is_indexed,
              'type': 'healpix'
             }, 
             upsert = True)


    def sphere2d_meta(self, sphere2d_key, is_indexed, pos_format = "geoJSON"):
        """
            Add metadata identifying the key with the source coordinates in 
            geoJSON/legacy pair formal. They will look like:
            
                { 
                    'key' : sphere2d_key
                    'is_indexed': is_indexed
                    'pos_format': pos_format
                }
            
            Paramters:
            ----------
                
                sphere2d_key: `str`
                    name of document keys identifying the position of the source in 
                    geoJSON or legacy pair format.
                
                is_indexed: `bool`
                    weather or not the srcs collection has been indexed using this key.
                
                pos_format: `str` 
                    format of the sphere2d_key, either geoJSON or legacy.
        """
        
        self.logger.info("adding metadata for 2dsphere coordinates named: %s"%sphere2d_key)
        self.db['meta'].update(
            {'_id': sphere2d_key},
            { 
              '_id' : sphere2d_key,
              'key' : sphere2d_key,
              'is_indexed': is_indexed,
              'pos_format': pos_format,
              'type': 'sphere2d'
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

