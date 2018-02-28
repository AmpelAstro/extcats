
extcats: tools to organize and query astronomical catalogs
==================================================================

This modules provides classes to import astronomical catalogs into 
a mongo noSql database, and to query this database for positional 
matches.

The preferred way to install the package is to download the tar 
from this Github page, unpack it in your favourite place, and 
run python setup.py install. This will also download some test
data to run the example notebooks. 

Alternatively, the pure python pacakge (witout examples) should be 
installable via pip.

The two main classes are:

    - CatalogPusher: will process the raw files with the catalog sources
    and creates a database.
    
    - CatalogQuery: will perform queries on the catalogs.

See examples and explanations in the notebook folder.
