import os
from sam2lca.acc2tax import get_mapping
from sam2lca.utils import get_json
from sam2lca.rocksdb import OPTS_read
import rocksdb
from shutil import rmtree
import pytest

test_dir = os.path.join(os.path.dirname(__file__))
db_dir = os.path.join(test_dir, "test_acc2tax_dbdir")
os.makedirs(db_dir, exist_ok=True)
maptype = "test"


def test_create_db():
    config = get_json("https://raw.githubusercontent.com/maxibor/sam2lca/master/data/acc2tax.json")
    get_mapping(map_config=config, maptype=maptype, dbdir=db_dir)


def test_acc2tax_db_query():
    config = get_json("https://raw.githubusercontent.com/maxibor/sam2lca/master/data/acc2tax.json")
    db = rocksdb.DB(
        os.path.join(db_dir, config["map_db"][maptype]),
        OPTS_read,
        read_only=True,
    )
    acc = "NC_000913.3"
    ret = db.get(acc.encode())
    assert ret.decode() == "511145"


def test_cleanup():
    rmtree(db_dir)
