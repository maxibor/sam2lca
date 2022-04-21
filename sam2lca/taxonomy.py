from pkg_resources import ExtractionError
from taxopy import TaxDb
from taxopy.exceptions import DownloadError
import pickle
from logging import error as logging_error
from os.path import isfile as os_isfile
from sys import exit as sys_exit
from pathlib import Path
from sam2lca.data.gtdb_taxdump import gtdb_taxdump
from sam2lca.utils import download_and_checksum, decompress
import os


def get_gtdb_taxonomy(dbdir, release="r207"):
    """Get GTDB taxonomy database
    Args:
        dbdir (str): Path to sam2lca database directory
        release (str, optional): GTDB release. Defaults to "r207".
    """

    gtdb_taxdump_local = dict()

    if release not in gtdb_taxdump:
        logging_error(f"GTDB release {release} not available")
        logging_error(f"Available releases: {', '.join(gtdb_taxdump.keys())}")
        sys_exit(1)
    else:
        for file in gtdb_taxdump[release]["url"]:
            fname = download_and_checksum(
                file_url=gtdb_taxdump[release]["url"][file],
                md5_url=gtdb_taxdump[release]["md5"][file],
                dbdir=dbdir,
            )
            gtdb_taxdump_local[file] = decompress(os.path.join(dbdir, fname))

    return gtdb_taxdump_local


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

    if db_type.lower() == "gtdb":
        gtdb_taxo = get_gtdb_taxonomy(dbdir)
        nodes = gtdb_taxo["nodes"]
        names = gtdb_taxo["names"]
        merged = gtdb_taxo["merged"]

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
