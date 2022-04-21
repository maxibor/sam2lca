import os
from sam2lca.data import acc2tax_default
from sam2lca.acc2tax import get_mapping
from sam2lca.rocksdb import OPTS_read
import rocksdb
from shutil import rmtree
import pytest

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
db_dir = os.path.join(test_dir, "test_acc2tax_dbdir")
os.makedirs(db_dir, exist_ok=True)
maptype = "test"


def test_create_db():
    get_mapping(map_config=acc2tax_default, maptype=maptype, dbdir=db_dir)


def test_acc2tax_db_query():
    db = rocksdb.DB(
        os.path.join(db_dir, acc2tax_default["map_db"][maptype]),
        OPTS_read,
        read_only=True,
    )
    acc = "NC_000913.3"
    ret = db.get(acc.encode())
    assert ret.decode() == "511145"


def test_cleanup():
    rmtree(db_dir)
