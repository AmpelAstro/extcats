
extcats: tools to organize and query astronomical catalogs
==================================================================

This modules provides classes to import astronomical catalogs into 
a mongo database, and to query this database for positional matches.

The three main classes are:

    - CatalogPusher: will process the raw files with the catalog sources
    and creates a database.
    
    - CatalogManager: will make the catalogs available to the query routines.
    
    - CatalogQuery: will perform queries on the catalogs.

See examples and explanations in the notebook folder.