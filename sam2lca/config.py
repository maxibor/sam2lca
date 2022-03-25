from taxopy import TaxDb
import rocksdb


def setup_taxopy_db(db_path, nodes=None, names=None, merged=None):
    """Setup taxopy database

    Args:
        db_path (str): Path to taxopy database
        nodes (str, optional): Path to nodes.dmp . Defaults to None.
        name (str, optional): Path to names.dmp. Defaults to None.
        merged (str, optional): Path to merged.dmp . Defaults to None.

    Returns:
        taxopy.TaxDb: Taxonomy database
    """
    return TaxDb(taxdb_dir=db_path, keep_files=True)


OPTS_create = rocksdb.Options()
OPTS_create.create_if_missing = True
OPTS_create.max_open_files = 250
OPTS_create.table_factory = rocksdb.BlockBasedTableFactory(
    filter_policy=rocksdb.BloomFilterPolicy(10),
    block_cache=rocksdb.LRUCache(6 * (1024**3)),
    block_cache_compressed=rocksdb.LRUCache(1.5 * (1024**3)),
)

OPTS_read = rocksdb.Options()
OPTS_read.create_if_missing = False
OPTS_read.max_open_files = 250
OPTS_read.table_factory = rocksdb.BlockBasedTableFactory(
    filter_policy=rocksdb.BloomFilterPolicy(10),
    block_cache=rocksdb.LRUCache(6 * (1024**3)),
    block_cache_compressed=rocksdb.LRUCache(1.5 * (1024**3)),
)
