import os
import pickle
import pytest
from shutil import rmtree
from sam2lca.main import update_database
from sam2lca.utils import taxid_to_lineage, count_reads_taxid

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
db_dir = os.path.join(test_dir, "test_utils_dbdir")
os.makedirs(db_dir, exist_ok=True)


read_ancestor_dict = {
    "read_A": 562,
    "read_B": 562,
    "read_C": 562,
    "read_D": 300267,
    "read_E": 300267,
    "read_F": 300267,
    "read_G": 300267,
    "read_H": 543,
    "read_I": 543,
    "read_J": 12908,  # unclassified taxid
}


@pytest.fixture(autouse=True)
def build_taxonomy_acc2tax_db(script_runner):
    update_database(
        dbdir=db_dir,
        acc2tax=None,
        taxonomy="test"
    )


@pytest.fixture(autouse=True)
def taxo_db():
    with open(os.path.join(db_dir, "test.pkl"), "rb") as tdb:
        taxo_db = pickle.load(tdb)
    return taxo_db


def test_taxid_to_lineage(taxo_db):
    taxid_counts = count_reads_taxid(read_ancestor_dict)
    test_result_dir = os.path.join(test_dir, "tmp")
    os.makedirs(test_result_dir, exist_ok=True)
    output = os.path.join(test_result_dir, "test.sam2lca")
    test_res = taxid_to_lineage(
        taxid_counts, output=output, process=2, nb_steps=7, taxo_db=taxo_db
    )
    assert test_res[12908]["count_taxon"] == 1
    assert test_res[543]["count_taxon"] == 2
    assert test_res[543]["count_descendant"] == 9
    assert test_res[1]["count_taxon"] == 0
    assert test_res[1]["count_descendant"] == 9
    assert test_res[300267]["count_taxon"] == 4
    assert test_res[300267]["count_descendant"] == 4


def test_cleanup():
    rmtree(db_dir)
    test_result_dir = os.path.join(test_dir, "tmp")
    rmtree(test_result_dir)
