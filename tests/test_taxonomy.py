import os
import pickle
from taxopy import Taxon
import pytest
from shutil import rmtree
from sam2lca.main import update_database
from sam2lca.lca import compute_lca

test_dir = os.path.join(os.path.dirname(__file__))
db_dir = os.path.join(test_dir, "test_taxonomy_dbdir")
os.makedirs(db_dir, exist_ok=True)


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


def test_taxo_db_query(taxo_db):
    taxon = Taxon(562, taxo_db)
    assert taxon.name == "Escherichia coli"


def test_compute_lca(taxo_db):
    read_dict = {
        "read_A": set({Taxon(562, taxo_db)}),
        "read_B": set({Taxon(562, taxo_db), Taxon(300267, taxo_db)}),
    }
    ancestors_dict = compute_lca(
        read_dict=read_dict, process=2, nb_steps=7, taxo_db=taxo_db
    )
    assert ancestors_dict["read_A"] == 562
    assert ancestors_dict["read_B"] == 543


def test_cleanup():
    rmtree(db_dir)
