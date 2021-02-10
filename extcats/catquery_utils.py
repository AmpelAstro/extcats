#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# collection of low-level functions used to query astronomical catalogs 
# for entries around a given target position. The catalogs are supposed to 
# be stored in mongodb collection.
#
# Author: M. Giomi (matteo.giomi@desy.de)


from numpy import append, unique, argmin
from math import radians, sqrt    # when running on scalars, math is much faster than numpy
from bson.code import Code
from healpy import ang2pix, pix2ang, get_all_neighbours, nside2resol
from astropy.table import Table
from astropy.coordinates import SkyCoord

import logging
module_logger = logging.getLogger(__name__)
logging.basicConfig(level = logging.INFO)

# javascript function to compute the angular distances
angdist_code = """
    function get_nearby(lat1=%s, lon1=%s, lat2=%.10f, lon2=%.10f) {
        var lat1_r = lat1 * Math.PI / 180.;
        var lat2_r = lat2 * Math.PI / 180.;
        var dLon_r = (lon2 - lon1) * Math.PI / 180.;
        var a = ( Math.sin(lat1_r)*Math.sin(lat2_r) + Math.cos(lat1_r)*Math.cos(lat2_r) * Math.cos(dLon_r));
        return (Math.acos(a) <= %.10f) ;
        }"""


def filters_logical_and(f1, f2, f3):
    """
        given three filters, returns the logical and of these filters
        in the order: f1 AND f2 AND f3.
        
        Parameters:
        -----------
            
            f[1/2/3]: `dict`
                filter expression to be combined.
        
        Returns:
        --------
            
            `dict` with the order-preserving combination of the filters.
    """
    
    if f1 is None and f3 is None:
        return f2
    elif f3 is None:
        return {'$and': [f1, f2]}
    elif f1 is None:
        return {'$and': [f2, f3]}
    else:
        return {'$and': [f1, f2, f3]}


def query_and_return(qfilter, projection, coll, find_one, to_table = True):
    """
        query a collection with a given filter and return the results.
        
        Parameters:
        -----------
        
            qfilter: `dict`
                query filter.
            
            coll: `pymongo.collection.Collection`
                collection to be queried.
            
            find_one: `bool`
                weather you want to use coll.find_one or coll.find.
            
            to_table: `bool`
                weather you want to return a list of dictionaries or an astropy Table.
        Returns:
        --------
            
            `list` or `astropy.table.Table` with query results, or None if no match.
    """
    
    if find_one:
        qresult = coll.find_one(qfilter, projection)
        if qresult is None:
            return None
        qresults = [qresult]
    else:
        qresults = list(coll.find(qfilter, projection))
    if len(qresults) == 0:
        return None
    elif to_table:
        return Table(qresults)
    else:
        return qresults


def searcharound_HEALPix(
    ra, dec, rs_arcsec, src_coll, hp_key, hp_order, hp_nest, hp_resol, 
    circular, ra_key, dec_key, find_one, sphere2d_key = None, pre_filter = None, post_filter = None, projection={"_id": 0}, logger  = None):
    """
        Returns sources in catalog contained in the the group of healpixels
        around the target coordinate that covers the search radius.
        
        NOTE that if circular is False, this function returns matches 
        in square-like pattern of HEALPixels that covers the search radius.
        
        Parameters:
        -----------
            
            ra/dec: `float`
                sky coordinates of the target direction. They have to be in
                the same format as those used to find the healpix id.
            
            rs_arcsec: `float`
                maximum allowed distance to target position (in arcsec).
            
            src_coll: `pymongo.collection.Collection`
                collection with the sources to be queried.
            
            hp_key: `str`
                name of the document key conatining the HEALPix index.
            
            hp_order: `int`
                order of the HEALPix grid used to group sources.
            
            hp_nest: `bool`
                whether or not the HEALPix grid has nested geometry.
            
            hp_resol: `float`
                resolution of the HEALPix map in arcsec, given by:
                    nside2resol(2**hp_order, arcmin = True) * 60.
            
            circular: `bool`, default: False
                if True, the results are skimmed removing sources outside of the
                circular search radius. In this case, the ra/dec keys arguments 
                needs to be specified.
            
            ra[dec]_key: `str`
                names of the coordinates in the catalog that will be used 
                for the queries.
            
            find_one: `bool`
                if True the collection is searched with the find_one method returning
                just the first result of the query. This behaviour is disabled if
                the circular flag is set (the HP query returns stuff within circle).

            pre_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain field are of interest:
                    pre_filter = {'field': field_number}
                This condition is checked BEFORE searching for healpix

            post_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain rcid are of interest:
                    post_filter = {'rcid': rcid_number}
                This condition is checked AFTER searching for healpix       
            
            logger: logging.logger istance
                logging for this function. If None, the module-level logger is used.
            
        Returns:
        --------
            
            cptable: `astropy.table.Table`/None
                astropy table of the catalog entry for the found counterpart. If
                no counpterpart is found returns None.
    """
    
    # start with the 9 neighbouring pixels on and around the target.
    nside = 2**hp_order
    target_pix = ang2pix(nside, ra, dec, nest = hp_nest, lonlat = True)
    neighbs = get_all_neighbours(nside, ra, dec, nest = hp_nest, lonlat = True)
    pix_group = append(target_pix, neighbs)
    pixgroup_scale = hp_resol*sqrt(len(pix_group))
    
    # increase the pixel group untill you cover the search radius.
    while pixgroup_scale/2. < rs_arcsec:
        pix_ras, pix_decs = pix2ang(nside, pix_group, lonlat = True, nest = hp_nest)
        new_pixels = get_all_neighbours(nside, pix_ras, pix_decs, nest = hp_nest, lonlat = True)
        pix_group = unique( append(pix_group, new_pixels.flatten() ) )
        pixgroup_scale = hp_resol*sqrt(len(pix_group))
    
    # cast a warning in case you have to enlarge the group too much
    if len(pix_group)>5000:
        if logger is None:
            logger = module_logger
        logger.warning(
            "search radius %.5f arcsec too big for healpixels of order %d (each ~ %.5f arcsec)"%
            (rs_arcsec, hp_order, hp_resol))
        logger.warning(
            "needed a group of %d healpixels with side of %d to cover search area."%
            (len(pix_group), sqrt(len(pix_group))))
    
    # remove non-existing neigbours (in case of E/W/N/S) and cast to list
    pix_group = pix_group[pix_group != -1].astype(int).tolist()
    
    # define your filter
    hp_query = {hp_key:{'$in': pix_group}}
    combined_filter = filters_logical_and(pre_filter, hp_query, post_filter)
    
    # query the database for sources in these pixels
    qresults = query_and_return(combined_filter, projection, src_coll, find_one, to_table=False)

    if qresults is None:
        return None
    if not circular:
        return Table(qresults)
    else:
        tab = Table(qresults)
        cols = set(tab.keys())
        if (ra_key not in cols or dec_key not in cols) and sphere2d_key is not None:
            ra_key, dec_key = "_ra", "_dec"
            tab[ra_key], tab[dec_key] = zip(*(q[sphere2d_key]["coordinates"] for q in qresults))
        dists = get_distances(ra = ra, dec = dec, table = tab, ra_key = ra_key, dec_key = dec_key)
        circular_tab = tab[dists <= rs_arcsec]
        if len(circular_tab) == 0:
            return None
        else:
            return circular_tab


def searcharound_9HEALPix(ra, dec, src_coll, hp_key, hp_order, hp_nest, hp_resol, 
    find_one, pre_filter = None, post_filter = None, projection={"_id": 0}):
    """
        Returns sources in catalog contained in the 9 healpixels
        around the target coordinate.
        
        Parameters:
        -----------
            
            ra/dec: `float`
                sky coordinates of the target direction. They have to be in
                the same format as those used to find the healpix id.
            
            src_coll: `pymongo.collection.Collection`
                collection with the sources to be queried.
            
            hp_key: `str`
                name of the document key conatining the HEALPix index.
            
            hp_order: `int`
                order of the HEALPix grid used to group sources.
            
            hp_nest: `bool`
                weather or not the HEALPix grid has nested geometry.
            
            find_one: `bool`
                if True the collection is searched with the find_one method returning
                just the first result of the query. if False, the method
                find is used, returning all matching documents. 

            pre_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain field are of interest:
                    pre_filter = {'field': field_number}
                This condition is checked BEFORE searching for healpix

            post_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain rcid are of interest:
                    post_filter = {'rcid': rcid_number}
                This condition is checked AFTER searching for healpix     
            
        Returns:
        --------
            
            cptable: `astropy.table.Table`/None
                astropy table of the catalog entry for the found counterpart. If
                no counpterpart is found returns None.
    """
    
    # find the index of the target pixel and its neighbours
    nside = 2**hp_order
    target_pix = int( ang2pix(nside, ra, dec, nest = hp_nest, lonlat = True) )
    neighbs = get_all_neighbours(nside, ra, dec, nest = hp_nest, lonlat = True)
    
    # remove non-existing neigbours (in case of E/W/N/S) and add center pixel
    pix_group = [int(pix_id) for pix_id in neighbs if pix_id != -1] + [target_pix]
    hp_query = {hp_key:{'$in': pix_group}}
    combined_filter = filters_logical_and(pre_filter, hp_query, post_filter)

    # query the database for sources in these pixels
    return query_and_return(combined_filter, projection, src_coll, find_one)


def searcharound_2Dsphere(ra, dec, rs_arcsec, src_coll, s2d_key, find_one, 
    pre_filter = None, post_filter = None, projection={"_id": 0}, logger = None):
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
            
            src_coll: `pymongo.collection.Collection`
                collection with the sources to be queried.
            
            s2d_key: `str`
                name of document key for the geoJSON/'legacy-pair' field.
            
            find_one: `bool`
                if True the collection is searched with the find_one method returning
                just the first result of the query. if False, the method
                find is used, returning all matching documents.

            pre_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain field are of interest:
                    pre_filter = {'field': field_number}
                This condition is checked BEFORE searching in the 2dsphere 

            post_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain rcid are of interest:
                    post_filter = {'rcid': rcid_number}
                This condition is checked AFTER searching in the 2dsphere     
            
            logger: logging.logger istance
                logging for this function. If None, the module-level logger is used.
            
        Returns:
        --------
            
            cptable: `astropy.table.Table`/None
                astropy table of the catalog entry for the found counterpart. If
                no counpterpart is found returns None.
    """
    
    # fold the RA between -180 and 180.
    if ra > 180:
        ra = ra - 360.
    
    # define filter
    geowithin = {"$geoWithin": { "$centerSphere": [[ra, dec], radians(rs_arcsec/3600.)]}}
    geoquery = {s2d_key: geowithin}
    combined_filter = filters_logical_and(pre_filter, geoquery, post_filter)

    # query and return
    tab = query_and_return(combined_filter, projection, src_coll, find_one)
    if tab is not None:
        ra_key, dec_key = "_ra", "_dec"
        tab[ra_key], tab[dec_key] = zip(*(q[s2d_key]["coordinates"] for q in tab))
    return tab


def searcharound_RAW(ra, dec, rs_arcsec, src_coll, ra_key, dec_key, find_one, box_scale = 2, pre_filter = None, post_filter = None, projection={"_id": 0, "pos": 0}):
    """
        Returns sources in catalog within rs_arcsec from target position.
        
        It first selects points within a box of radius box_scale times larger than the
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
            
            src_coll: `pymongo.collection.Collection`
                collection with the sources to be queried.
            
            ra[dec]_key: `str`
                names of the coordinates in the catalog that will be used 
                for the queries.
            
            box_scale: `float`
                size of the square used to pre-search candidates, in units of rs_arcsec.
            
            find_one: `bool`
                if True the collection is searched with the find_one method returning
                just the first result of the query. if False, the method
                find is used, returning all matching documents. 

            pre_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain field are of interest:
                    pre_filter = {'field': field_number}
                This condition is checked BEFORE searching in the box 

            post_filter: 'dict'
                contains an additional condition on the queried objects.
                For example, if only stars in a certain rcid are of interest:
                    post_filter = {'rcid': rcid_number}
                This condition is checked AFTER searching in the box 
            
        Returns:
        --------
            
            cptable: `astropy.table.Table`/None
                astropy table of the catalog entry for the found counterpart. If
                no counpterpart is found returns None.
    """
    
    # wrap your angular distance query into into javascript 
    query_func = Code(angdist_code%(
        "this.%s"%dec_key, "this.%s"%ra_key, dec, ra, radians(rs_arcsec/3600.)))
    
    # restrict first to a box in the query
    box_size = box_scale * (rs_arcsec/3600.)
    qfilter =  { 
            ra_key : { "$gte" :  ra-box_size, "$lte" : ra+box_size}, 
            dec_key: { "$gte" :  dec-box_size, "$lte" : dec+box_size},
            "$where": query_func
            }
    combined_filter = filters_logical_and(pre_filter, qfilter, post_filter)

    # query and return
    return query_and_return(combined_filter, projection, src_coll, find_one)


def get_distances(ra, dec, table, ra_key, dec_key):
    """
        given an astropy Table containing sources with position keys ra[dec]_key, 
        return an array with the angluar distance (in arcsec) of each table entry
        to the target ra, dec coordinates.
        
        Parameters:
        -----------
            
            ra/dec: `float`
                sky coordinates (deg), of the target direction. They have to be in the 
                same refernce system as the ra[dec]_key in the Table.
            
            table: `astropy.table.Table`
                table with the sources among which the closest to target will be
                selected.
            
            ra[dec]_key: `str`
                names of the coordinates (in deg) in the input table.
        
        Returns:
        --------
            
            dist: `array` 
                angular distances (in arcsec) each table entry to the target.
    """
    
    target = SkyCoord(ra, dec, unit = "deg")
    matches_pos = SkyCoord(table[ra_key], table[dec_key], unit = "deg")
    return target.separation(matches_pos).arcsecond


def get_closest(ra, dec, table, ra_key, dec_key):
    """
        given an astropy Table containing sources with position keys ra[dec]_key, 
        return the entry closest to the target.
        
        Parameters:
        -----------
            
            ra/dec: `float`
                sky coordinates (J2000, in degrees) of the target direction.
            
            table: `astropy.table.Table`
                table with the sources among which the closest to target will be
                selected.
            
            ra[dec]_key: `str`
                names of the coordinates in the input table.
        
        Returns:
        --------
            
            cptable: `astropy.table.Table`
                closest source found in table.
            
            dist: `float`/None
                angular distance (in arcsec) of the closest source to the target.
    """
    d2t = get_distances(ra, dec, table, ra_key, dec_key)
    match_id = argmin(d2t)
    return table[match_id], d2t[match_id]


def random_point_sphere(npoints, rnd_seed = None):
    """
        Utitlity function to provide test coordinates to measure 
        query times. It computes randomly distributed points on 
        a 2d sphere.
        
        Parameters:
        -----------
        
        npoints: `int`
            number of points to generate.
        
        rnd_seed: `float` or None
            if not None, this seed will be passed to np.random. Else use numpy 
            defaults.
        
        Returns:
        --------
        
        numpy array structured as [[ra1, dec1], [ra2, dec2], ...]
    """
    import numpy as np
    if not rnd_seed is None:
        np.random.seed(rnd_seed)
    points = np.zeros(2*npoints).reshape(npoints, 2)
    phi = - np.pi + 2*np.pi*np.random.ranf(npoints)
    theta = np.arcsin( 2*np.random.ranf(npoints) - 1 )
    points = np.stack( (np.degrees(phi), np.degrees(theta)), axis = 1)
    return points
    

