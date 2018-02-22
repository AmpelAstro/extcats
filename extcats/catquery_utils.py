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

from extcats import CatalogManager as CM

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level = logging.INFO)

def searcharound_HEALPix(ra, dec, rs_arcsec, src_coll, hp_key, hp_order, hp_nest, hp_resol):
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
            
            src_coll: `pymongo.collection.Collection`
                collection with the sources to be queried.
            
            hp_key: `str`
                name of the document key conatining the HEALPix index.
            
            hp_order: `int`
                order of the HEALPix grid used to group sources.
            
            hp_nest: `bool`
                weather or not the HEALPix grid has nested geometry.
            
            hp_resol: `float`
                resolution of the HEALPix map in arcsec, given by:
                    nside2resol(2**hp_order, arcmin = True) * 60.
                
            
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
    
    # cast a warning in case you have to enlarge the group
    if pixgroup_scale/2. < rs_arcsec:
        logging.warning(
        "search radius %.5f arcsec too big for healpixels of order %d (each ~ %.5f arcsec)"%
        (rs_arcsec, hp_order, hp_resol))
    
    # increase the pixel group untill you cover the search radius.
    while pixgroup_scale/2. < rs_arcsec:
        pix_ras, pix_decs = pix2ang(nside, pix_group, lonlat = True, nest = hp_nest)
        new_pixels = get_all_neighbours(nside, pix_ras, pix_decs, nest = hp_nest, lonlat = True)
        pix_group = unique( append(pix_group, new_pixels.flatten() ) )
        pixgroup_scale = hp_resol*sqrt(len(pix_group))
        big_group = True
    if len(pix_group)>9:
        logging.warning(
        "needed a group of %d healpixels with side of %d to cover search area."%(len(pix_group), sqrt(len(pix_group))))
    
    # remove non-existing neigbours (in case of E/W/N/S) and cast to list
    pix_group = pix_group[pix_group != -1].astype(int).tolist()
    
    # query the database for sources in these pixels
    qresults = [o for o in src_coll.find( {hp_key: { "$in": pix_group }} )]
    if len(qresults)==0:
        return None
    else:
        return Table(qresults)


def searcharound_9HEALPix(ra, dec, src_coll, hp_key, hp_order, hp_nest, hp_resol):
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
    
    # query the database for sources in these pixels
    qresults = [o for o in src_coll.find( {hp_key: { "$in": pix_group }} )]
    if len(qresults)==0:
        return None
    else:
        return Table(qresults)


def searcharound_2Dsphere(ra, dec, rs_arcsec, src_coll, s2d_key):
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
            
        Returns:
        --------
            
            cptable: `astropy.table.Table`/None
                astropy table of the catalog entry for the found counterpart. If
                no counpterpart is found returns None.
    """
    
    # fold the RA between -180 and 180.
    if ra > 180:
        logging.warning(
        "2dsphere queries needs R.A in [-180, 180] deg. Folding with: ra = ra - 360 if ra > 180.")
        ra = ra - 360.
    
    # query and return
    geowithin={"$geoWithin": { "$centerSphere": [[ra, dec], radians(rs_arcsec/3600.)]}}
    qresults = [o for o in src_coll.find({s2d_key: geowithin})]
    if len(qresults)==0:
        return None
    else:
        return Table(qresults)


def searcharound_RAW(ra, dec, rs_arcsec, src_coll, ra_key, dec_key, box_scale = 2):
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
            
        Returns:
        --------
            
            cptable: `astropy.table.Table`/None
                astropy table of the catalog entry for the found counterpart. If
                no counpterpart is found returns None.
    """
    
    # wrap your angular distance query into into javascript 
    angdist_code = """
    function get_nearby(lat1=%s, lon1=%s, lat2=%.10f, lon2=%.10f) {
        var lat1_r = lat1 * Math.PI / 180.;
        var lat2_r = lat2 * Math.PI / 180.;
        var dLon_r = (lon2 - lon1) * Math.PI / 180.;
        var a = ( Math.sin(lat1_r)*Math.sin(lat2_r) + Math.cos(lat1_r)*Math.cos(lat2_r) * Math.cos(dLon_r));
        return (Math.acos(a) <= %.10f) ;
        }"""
    query_func = Code(angdist_code%(
        "this.%s"%dec_key, "this.%s"%ra_key, dec, ra, radians(rs_arcsec/3600.)))
    
    # restrict first to a box in the query
    box_size = box_scale * (rs_arcsec/3600.)
    qfilter =  { 
            ra_key : { "$gte" :  ra-box_size, "$lte" : ra+box_size}, 
            dec_key: { "$gte" :  dec-box_size, "$lte" : dec+box_size},
            "$where": query_func
            }
    qresults = [o for o in src_coll.find(qfilter)]
    if len(qresults)==0:
        return None
    else:
        return Table(qresults)


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
    
    target = SkyCoord(ra, dec, unit = "deg")
    matches_pos = SkyCoord(table[ra_key], table[dec_key], unit = "deg")
    d2t = target.separation(matches_pos).arcsecond
    match_id = argmin(d2t)
    return table[match_id], d2t[match_id]


def random_point_sphere(npoints):
    """
        Utitlity function to provide test coordinates to measure 
        query times. It computes randomly distributed points on 
        a 2d sphere.
        
        Parameters:
        -----------
        
        npoints: `int`
            number of points to generate.
        
        Returns:
        --------
        
        numpy array structured as [[ra1, dec1], [ra2, dec2], ...]
    """
    import numpy as np
    np.random.seed(42)
    points = np.zeros(2*npoints).reshape(npoints, 2)
    phi = - np.pi + 2*np.pi*np.random.ranf(npoints)
    theta = np.arcsin( 2*np.random.ranf(npoints) - 1 )
    points = np.stack( (np.degrees(phi), np.degrees(theta)), axis = 1)
    return points
    

