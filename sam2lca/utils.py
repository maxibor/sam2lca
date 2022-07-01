import os
import json
import rocksdb
from taxopy import Taxon
import pandas as pd
import logging
import warnings
from collections import Counter
from tqdm.contrib.concurrent import thread_map
from tqdm import tqdm
from functools import partial
import hashlib
import urllib.request as urllib
from urllib.error import URLError
import shutil
from subprocess import check_output
import sys
import gzip

warnings.simplefilter(action="ignore", category=FutureWarning)

taxo_ranks = [
    'strain',
    'species',
    'genus',
    'family',
    'order',
    'class',
    'phylum',
    'superkingdom'
]

def count_reads_taxid(read_taxid_dict):
    """Returns number of reads matching TAXID

    Args:
        read_taxid_dict (dict): {read_name(str): TAXID(int)}
    Returns:
        {TAXID(int): count(int)}
    """
    return dict(Counter(read_taxid_dict.values()))


def output_file(thepath):
    out = os.path.basename(thepath).split(".")[:-1]
    if len(out) == 0:
        out = {"sam2lca": f"{thepath}.sam2lca", "bam": f"{thepath}.bam"}
    else:
        out = {
            "sam2lca": ".".join(out) + ".sam2lca",
            "bam": ".".join(out) + ".sam2lca.bam",
        }

    return out


def check_extension(filename):
    """Check alignment file format to give correct open mode

    Args:
        filename (str): Path to alignment file

    Returns:
        str: opening mode

    Raises:
        Exception: Extension not supported
    """
    extension = filename.split(".")[-1]
    modes = {"bam": "rb", "sam": "r", "cram": "rc"}
    try:
        return modes[extension]
    except KeyError:
        raise Exception(f"{extension} file extension not supported")


def reassign_count_lineage(taxid_items, read_taxid_dict_reassign, taxo_db):
    """Compute total counts of reads matching to taxons and their descendants

    Args:
        taxid_items (dict_item): (TAXID(int), read_count(int))
    """
    try:
        taxid, read_count = taxid_items
        taxon = Taxon(taxid, taxo_db)
        lineage = taxon.taxid_lineage
    except Exception:
        read_taxid_dict_reassign[taxid] = [read_count, read_count]
        return
    if taxid not in read_taxid_dict_reassign:
        read_taxid_dict_reassign[taxid] = [read_count, 0]
    else:
        read_taxid_dict_reassign[taxid][0] = read_count
    for t in lineage:
        if t not in read_taxid_dict_reassign:
            read_taxid_dict_reassign[t] = [0, read_count]
        else:
            read_taxid_dict_reassign[t][1] += read_count


def taxid_to_lineage_single(taxid_items, taxid_info_dict, taxo_db):
    """Retrieve Taxonomic lineage and species name from NCBI TAXID

    Args:
        taxid_items (dict_item): (TAXID(int), read_count(int))
    Returns:
        {taxid(int): {"name":taxon_name(str), "rank":taxon_rank(str), "count":read_count(int), "lineage":taxon_lineage(list)}}
    """
    try:
        taxid, read_counts = taxid_items
        read_count_taxon, read_count_descendant = read_counts
        taxon = Taxon(taxid, taxo_db)
        sciname = taxon.name
        rank = taxon.rank
        lineage = taxon.rank_name_dictionary
        taxid_lineage = taxon.taxid_lineage
    except Exception:
        sciname = "unclassified sequences"
        rank = "NA"
        lineage = "NA"
        taxid_lineage = "NA"

    taxid_info_dict[taxid] = {
        "name": sciname,
        "rank": rank,
        "count_taxon": read_count_taxon,
        "count_descendant": read_count_descendant,
        "lineage": lineage,
        "lineage_taxid": taxid_lineage,
    }


def taxid_to_lineage(taxid_count_dict, output, process, nb_steps, taxo_db):
    logging.info(
        f"Step {6 if nb_steps == 7 else 7 }/{nb_steps}: Converting TAXIDs to taxonomic lineages"
    )

    read_taxid_dict_reassign = dict()
    reassign_count_lineage_partial = partial(
        reassign_count_lineage,
        read_taxid_dict_reassign=read_taxid_dict_reassign,
        taxo_db=taxo_db,
    )
    thread_map(
        reassign_count_lineage_partial,
        taxid_count_dict.items(),
        chunksize=1,
        max_workers=process,
    )
    taxid_info_dict = dict()
    taxid_to_lineage_single_partial = partial(
        taxid_to_lineage_single, taxid_info_dict=taxid_info_dict, taxo_db=taxo_db
    )
    res = thread_map(
        taxid_to_lineage_single_partial,
        read_taxid_dict_reassign.items(),
        chunksize=1,
        max_workers=process,
    )

    df = pd.DataFrame(taxid_info_dict).transpose()
    df["lineage"] = (
        df["lineage"]
        .astype(str)
        .str.replace("[\[\]\{\}]", "")
        .str.replace(", ", " || ")
        .str.replace("'", "")
    )
    df.sort_values("count_descendant", inplace=True, ascending=False)
    df.to_csv(f"{output}.csv", index_label="TAXID")

    with open(f"{output}.json", "w") as write_file:
        json.dump(taxid_info_dict, write_file)
    logging.info(
        f"Step {7 if nb_steps == 7 else 8 }/{nb_steps}: writing sam2lca results:\n* JSON to {output}.json\n* CSV to {output}.csv"
    )

    return taxid_info_dict


def filter_taxid_info_dict(taxid_info_dict, rank, minread):
    taxid_info_dict_filtered = dict()
    for entry in taxid_info_dict:
        if taxid_info_dict[entry]['rank'] == rank and taxid_info_dict[entry]['count_descendant'] >= minread:
            taxid_info_dict_filtered[entry] = taxid_info_dict[entry]
    return taxid_info_dict_filtered

def get_db_connection(db_path):
    return rocksdb.DB(db_path, opts=rocksdb.Options(), read_only=True)


class TqdmUpTo(tqdm):
    """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""

    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)  # will also set self.n = b * bsize


def wc(filename, compressed=True):
    if compressed:
        cmd = f"zcat < {filename} | wc -l"
    else:
        cmd = f"wc -l {filename}"
    return int(check_output(cmd, shell=True).split()[0])


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def download_and_checksum(file_url, md5_url, dbdir):
    """Download  file

    Args:
        file_url (str): URL of  file
        md5_url (str): URL of md5  file
        dbdir (str): path to  local directory
    Returns:
        fname (str):  filename
    """

    md5_fname = os.path.basename(md5_url)
    fname = os.path.basename(file_url)
    if file_url.startswith("http") or file_url.startswith("ftp"):
        urllib.urlretrieve(md5_url, f"{os.path.join(dbdir, md5_fname)}")

        with TqdmUpTo(unit="B", unit_scale=True, miniters=1, desc=fname) as t:
            urllib.urlretrieve(
                file_url,
                filename=f"{os.path.join(dbdir, fname)}",
                reporthook=t.update_to,
                data=None,
            )
    else:
        shutil.copy(file_url, f"{os.path.join(dbdir, fname)}")
        shutil.copy(md5_url, f"{os.path.join(dbdir, md5_fname)}")
    with open(f"{os.path.join(dbdir,md5_fname)}", "r") as hf:
        for line in hf:
            md5_hash = line.rstrip().split()[0]
    file_hash = md5(f"{os.path.join(dbdir, fname)}")
    try:
        assert file_hash == md5_hash
    except AssertionError:
        logging.error(f"Error while downloading {file_url}, please try again")
        sys.exit(1)
    return fname


def decompress(filepath):
    basename = os.path.basename(filepath).replace(".gz", "")
    decomp_filepath = os.path.join(os.path.dirname(filepath), basename)
    with gzip.open(filepath, "rb") as f_in:
        with open(decomp_filepath, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decomp_filepath

def get_json(json_path):
    """Reads map_config.json file
    Args:
        json_path (str): Path to acc2tax.json file
    Returns:
        dict: json
    """
    
    if json_path.startswith("http"):
        try:
            with urllib.urlopen(json_path) as f:
                config = json.loads(f.read())
        except URLError as e:
            logging.error(e)
            sys.exit(1)
    else:
        with open(json_path, 'r') as f:
            config = json.load(f)        

    return config
