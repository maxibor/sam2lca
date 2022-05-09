import os
from shutil import rmtree
from sam2lca.main import list_available_db
import pytest

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
db_dir = os.path.join(test_dir, "test_update_dbdir")
os.makedirs(db_dir, exist_ok=True)


def test_help(script_runner):
    ret = script_runner.run("sam2lca", "update-db", "--help")
    assert ret.success


def test_build_taxonomy_acc2tax_db(script_runner):
    ret = script_runner.run(
        "sam2lca",
        "--dbdir",
        db_dir,
        "update-db",
        "--taxonomy",
        "test",
        "--acc2tax",
        "test",
    )
    assert ret.success


def test_list_db_help(script_runner):
    ret = script_runner.run("sam2lca", "list-db", "--help")
    assert ret.success


def test_list_db_cli(script_runner):
    ret = script_runner.run("sam2lca", "--dbdir", db_dir, "list-db")
    assert ret.success


def test_list_db_api():
    taxo, acc2tax = list_available_db(db_dir)
    assert taxo == ["test"]
    assert acc2tax == ["test"]


def test_cleanup():
    rmtree(db_dir)
