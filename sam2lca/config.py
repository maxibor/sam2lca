import ete3
import rocksdb

NCBI = ete3.NCBITaxa()

OPTS_create = rocksdb.Options()
OPTS_create.create_if_missing = True
OPTS_create.max_open_files = 250
OPTS_create.table_factory = rocksdb.BlockBasedTableFactory(
    filter_policy=rocksdb.BloomFilterPolicy(10),
    block_cache=rocksdb.LRUCache(6 * (1024 ** 3)),
    block_cache_compressed=rocksdb.LRUCache(1.5 * (1024 ** 3)))

OPTS_read = rocksdb.Options()
OPTS_read.create_if_missing = False
OPTS_read.max_open_files = 250
OPTS_read.table_factory = rocksdb.BlockBasedTableFactory(
    filter_policy=rocksdb.BloomFilterPolicy(10),
    block_cache=rocksdb.LRUCache(6 * (1024 ** 3)),
    block_cache_compressed=rocksdb.LRUCache(1.5 * (1024 ** 3)))