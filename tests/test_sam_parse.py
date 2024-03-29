import os
from sam2lca.parse_alignment import Alignment
from sam2lca.utils import get_json
from shutil import rmtree
import pickle
import pytest

test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, "data")
db_dir = os.path.join(test_dir, "test_sam_parse_dbdir")
os.makedirs(db_dir, exist_ok=True)
maptype = "test"




def test_build_taxonomy_acc2tax_db(script_runner):
    ret = script_runner.run(
        "sam2lca",
        "--dbdir",
        db_dir,
        "update-db",
        "--taxonomy",
        f"{maptype}",
        "--acc2tax",
        f"{maptype}",
    )
    assert ret.success

def test_parse_sam():

    with open(os.path.join(db_dir, f"{maptype}.pkl"), "rb") as tdb:
        taxo_db = pickle.load(tdb)

    sam = os.path.join(data_dir, "aligned.sorted.bam")
    al = Alignment(al_file=sam, nb_steps=7)
    config = get_json("https://raw.githubusercontent.com/maxibor/sam2lca/master/data/acc2tax.json")
    al.get_refs_taxid(acc2tax_dbname = os.path.join(db_dir, config["map_db"][maptype]), taxo_db=taxo_db)
    
    read_dict_i_0_8 = al.get_reads(
        process=2,
        identity=0.8,
        edit_distance=None,
        minlength=30,
        check_conserved=False,
    )
    assert len(read_dict_i_0_8["escherichia_coli_1057"]) == 2

    read_dict_i_1 = al.get_reads(
        process=2,
        identity=1,
        edit_distance=None,
        minlength=30,
        check_conserved=False,
    )
    assert len(read_dict_i_1["escherichia_coli_1057"]) == 1

    read_dict_d_2 = al.get_reads(
        process=2,
        identity=1,
        edit_distance=2,
        minlength=30,
        check_conserved=False,
    )
    assert len(read_dict_d_2["escherichia_coli_1057"]) == 2

    read_dict_d_0 = al.get_reads(
        process=2,
        identity=1,
        edit_distance=0,
        minlength=30,
        check_conserved=False,
    )
    assert len(read_dict_d_0["escherichia_coli_1057"]) == 1



def test_cleanup():
    rmtree(db_dir)