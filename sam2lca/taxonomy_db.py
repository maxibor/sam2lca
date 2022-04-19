from pkg_resources import ExtractionError
from taxopy import TaxDb
from taxopy.exceptions import DownloadError
import pickle
from logging import error as logging_error
from os.path import isfile as os_isfile
from sys import exit as sys_exit
from pathlib import Path


def setup_taxopy_db(dbdir, db_type="ncbi", nodes=None, names=None, merged=None):
    """Setup taxopy database

    Args:
        dbdir (str): Path to sam2lca database directory
        db_type(str): Type of taxonomy database, ncbi, or custom
        nodes (str, optional): Path to nodes.dmp . Defaults to None.
        name (str, optional): Path to names.dmp. Defaults to None.
        merged (str, optional): Path to merged.dmp . Defaults to None.

    Returns:
        taxopy.TaxDb: Taxonomy database
    """
    files_to_clean = [
        "merged.dmp",
        "names.dmp",
        f"{db_type.lower()}.pkl",
        "nodes.dmp",
        "taxdump.tar.gz",
    ]
    for file in files_to_clean:
        try:
            Path(f"{dbdir}/{file}").unlink()
        except FileNotFoundError:
            continue
    if db_type.lower() != "ncbi" and not nodes and not names and not merged:
        logging_error(
            "Custom taxonomy database cannot be created without specifying nodes.dmp, names.dmp, and merged.dmp files"
        )
        sys_exit(1)
    try:
        kf = True if db_type.lower() != "ncbi" else False
        TAXDB = TaxDb(
            taxdb_dir=dbdir,
            keep_files=kf,
            nodes_dmp=nodes,
            names_dmp=names,
            merged_dmp=merged,
        )
        Path(dbdir).mkdir(exist_ok=True)
        with open(f"{dbdir}/{db_type}.pkl", "wb") as f:
            pickle.dump(TAXDB, f)
    except (DownloadError, ExtractionError) as e:
        logging_error("Error while creating taxonomy database")
        logging_error(e)
        sys_exit(1)


def load_taxonomy_db(dbdir, db_type="ncbi"):
    """Load taxonomy database

    Args:
        dbdir (str): Path to sam2lca database directory
        db_type(str): Type of taxonomy database, NCBI, GTDB, or custom

    Returns:
        taxopy.TaxDb: Taxonomy database
    """
    db_type = db_type.lower()
    if os_isfile(f"{dbdir}/{db_type}.pkl"):
        with open(f"{dbdir}/{db_type}.pkl", "rb") as f:
            TAXDB = pickle.load(f)
        return TAXDB
    else:
        logging_error(
            f"Taxonomy database {db_type} not found, please run the update_db command"
        )
        sys_exit(1)
