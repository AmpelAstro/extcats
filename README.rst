*******
extcats
*******

tools to organize and query astronomical catalogs
#################################################


This modules provides classes to import astronomical catalogs into 
a **mongodb** database, and to efficiently query this database for 
positional matches.


Description:
############

The two main classes of this module are:

    - **CatalogPusher**: will process the raw files with the catalog sources and creates a database. See *insert_example* notebook for more details and usage instruction.
    
    - **CatalogQuery**: will perform queries on the catalogs. See *query_example* for examples and benchmarking.

Supported queries includes:

 - all the sources with a certain distance.
 - closest source at a given position.
 - binary search: return yes/no if anything is around the positon.
 - user defined queries.

The first item on the above list (cone search around target) provides the basic block for the other two types of positional-based queries. The code supports tree types of basic
cone-search queries, depending on the indexing strategy of the database.

    - using **HEALPix**: if the catalog sources have been assigned an HEALPix index (using `healpy <https://healpy.readthedocs.io/en/latest/#>`_).
     
    - using **GeoJSON** (or 'legacy coordinates'): if the catalog documents have the 
      position arranged in one of these two formats (`example 
      <https://docs.mongodb.com/manual/geospatial-queries/>`_), the query is based on
      the ``$geoWithin`` and ``$centerSphere`` mongo operators.
    
    - **raw**: this method uses the ``$where`` keyword to evaluate on each document a ``javascript``
      function computing the angular distance between each source and the target. This method 
      does not require any additional field to be added to the catalog but has, in general, 
      poorer performances with respect to the methods above.
      
All the core functions are defined in the ``catquery_utils`` module. In all cases the 
results of the queries will be return an ``astropy.table.Table`` objects.


Notes on indexing and query performances:
-----------------------------------------

The recommended method to index and query catalogs is based on the GeoJSON coorinate type.
See the *example_insert* notebook for how this can be implemented. 


Performant queries requires the database indexes to reside in the RAM. The indexes are 
efficiently compressed by mongodb default engine (WiredTiger), however there is little
redundant (and hence compressible) information in accurately measured coordinate pairs.
As a consequence, GeoJSON type indexes tends to require fair amount of free memory (of 
the order 40 MB for 2M entries). For large catalogs (and / or small RAM) indexing on 
coordinates might not be feasible. In this case, the HEALPix based indexing should 
be used. As (possibly) many sources shares the same HEALPix index, compression is 
more efficient into moderating RAM usage.

Installation:
^^^^^^^^^^^^^

The preferred way to install the package is to download the tar 
from this Github, unpack it to your favourite place, and then run:
::

    python setup.py install

This will also download some test data to run the example notebooks.

Alternatively, the pure python module (without notebooks examples
or test data) is installable via pip:
::

    pip install extcats


Usefull links:
--------------

 - `mongodb installation <https://docs.mongodb.com/manual/administration/install-community/>`_
 - `healpy <https://healpy.readthedocs.io/en/latest/#>`_
 - `astropy <http://www.astropy.org/>`_
