import os
import json
import pytest
from shutil import rmtree
from sam2lca import sam2lca

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
test_result_dir = os.path.join(test_dir, "tmp")
os.makedirs(test_result_dir, exist_ok=True)
db_dir = os.path.join(test_dir, "test_acc2tax_dbdir")
os.makedirs(db_dir, exist_ok=True)

output = os.path.join(test_dir, "tmp", "test")


def test_help(script_runner):
    ret = script_runner.run("sam2lca", "analyze", "--help")
    assert ret.success


@pytest.fixture(autouse=True)
def test_build_taxonomy_acc2tax_db(script_runner):
    nodes = os.path.join(data_dir, "taxonomy", "nodes.dmp")
    names = os.path.join(data_dir, "taxonomy", "names.dmp")
    merged = os.path.join(data_dir, "taxonomy", "merged.dmp")
    ret = script_runner.run(
        "sam2lca",
        "--dbdir",
        db_dir,
        "update-db",
        "--taxonomy",
        "test",
        "--taxo_nodes",
        nodes,
        "--taxo_names",
        names,
        "--taxo_merged",
        merged,
        "--acc2tax",
        "test",
    )


def test_analyze_cli(script_runner):
    bam = os.path.join(data_dir, "aligned.sorted.bam")
    ret = script_runner.run(
        "sam2lca",
        "--dbdir",
        db_dir,
        "analyze",
        "--taxonomy",
        "test",
        "--acc2tax",
        "test",
        "--output",
        output,
        bam,
    )
    assert ret.success


def test_analyze_results():
    json_file = os.path.join(data_dir, "results", "aligned.sorted.sam2lca.json")
    bam = os.path.join(data_dir, "aligned.sorted.bam")
    os.makedirs(test_result_dir, exist_ok=True)
    with open(json_file, "r") as f:
        target_results = json.load(f)
    sam2lca(sam=bam, dbdir=db_dir, taxonomy="test", acc2tax="test", output=output)
    with open(f"{output}.sam2lca.json", "r") as f:
        test_results = json.load(f)
    assert target_results == test_results

def test_analyze_distance_cli(script_runner):
    bam = os.path.join(data_dir, "aligned.sorted.bam")
    ret = script_runner.run(
        "sam2lca",
        "--dbdir",
        db_dir,
        "analyze",
        "--taxonomy",
        "test",
        "--acc2tax",
        "test",
        "--distance",
        "1",
        "--output",
        output,
        bam,
    )
    assert ret.success


def test_cleanup():
    rmtree(db_dir)
    rmtree(test_result_dir)
