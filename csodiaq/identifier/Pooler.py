from csodiaq.loaders.library.LibraryLoaderContext import LibraryLoaderContext
from csodiaq.loaders.query.QueryLoaderContext import QueryLoaderContext

class Pooler:
    def __init__(self, libraryFile, queryFile):
        self.libraryDict = LibraryLoaderContext(libraryFile).load_csodiaq_library_dict()
        self.queryContext = QueryLoaderContext(queryFile)