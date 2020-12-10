import ete3
import rocksdb

NCBI = ete3.NCBITaxa()

OPTS = rocksdb.Options()
OPTS.table_factory = rocksdb.BlockBasedTableFactory(
    filter_policy=rocksdb.BloomFilterPolicy(10),
    block_cache=rocksdb.LRUCache(6 * (1024 ** 3)),
    block_cache_compressed=rocksdb.LRUCache(1.5 * (1024 ** 3)))