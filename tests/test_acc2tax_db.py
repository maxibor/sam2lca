import os
from sam2lca.acc2tax import base_map_config
from sam2lca.acc2tax_rocksdb import get_mapping
from sam2lca.rocksdb_config import OPTS_read
import rocksdb
from shutil import rmtree
import pytest

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
db_dir = os.path.join(test_dir, "test_acc2tax_dbdir")
os.makedirs(db_dir, exist_ok=True)
maptype = "test"


def test_create_db():
    get_mapping(map_config=base_map_config, maptype=maptype, dbdir=db_dir)


def test_acc2tax_db_query():
    db = rocksdb.DB(
        os.path.join(db_dir, base_map_config["map_db"][maptype]),
        OPTS_read,
        read_only=True,
    )
    acc = "NC_000913.3"
    ret = db.get(acc.encode())
    assert ret.decode() == "511145"


def test_cleanup():
    rmtree(db_dir)
