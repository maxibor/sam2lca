import rocksdb

OPTS_create = rocksdb.Options()
OPTS_create.create_if_missing = True
OPTS_create.max_open_files = 500
OPTS_create.target_file_size_base = 10**8
OPTS_create.target_file_size_multiplier = 3

OPTS_read = rocksdb.Options()
OPTS_read.create_if_missing = False
OPTS_read.max_open_files = 500
OPTS_read.target_file_size_base = 10**8
OPTS_read.target_file_size_multiplier = 3
