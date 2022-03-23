import os
import json
import rocksdb
from sam2lca.config import NCBI
import pandas as pd
import logging
import warnings
from collections import Counter, ChainMap
from tqdm.contrib.concurrent import process_map

warnings.simplefilter(action="ignore", category=FutureWarning)


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


def taxid_to_lineage_single(taxid_items):
    """Retrieve Taxonomic lineage and species name from NCBI TAXID

    Args:
        taxid_items (dict_item): (TAXID(int), read_count(int))
    Returns:
        {taxid(int): {"name":taxon_name(str), "rank":taxon_rank(str), "count":read_count(int), "lineage":taxon_lineage(list)}}
    """
    try:
        taxid, read_count = taxid_items
        sciname = NCBI.get_taxid_translator([taxid])[taxid]
        rank = NCBI.get_rank([taxid])[taxid]
        taxid_lineage = NCBI.get_lineage(taxid)
        scinames = NCBI.get_taxid_translator(taxid_lineage)
        ranks = NCBI.get_rank(taxid_lineage)
        lineage = [{ranks[taxid]: scinames[taxid]} for taxid in taxid_lineage]
    except Exception:
        sciname = "NA"
        rank = "NA"
        lineage = "NA"
    return {
        taxid: {
            "name": sciname,
            "rank": rank,
            "count": read_count,
            "lineage": lineage,
        }
    }


def taxid_to_lineage(taxid_count_dict, output, process, nb_steps):
    logging.info(
        f"Step {6 if nb_steps == 7 else 7 }/{nb_steps}: Converting TAXIDs to taxonomic lineages"
    )
    res = process_map(
        taxid_to_lineage_single,
        taxid_count_dict.items(),
        chunksize=1,
        max_workers=process,
    )
    res = dict(ChainMap(*res))

    df = pd.DataFrame(res).transpose().sort_values("count", ascending=False)
    df["lineage"] = (
        df["lineage"].astype(str).str.replace("[\[\]\{\}]", "").str.replace(", ", " - ")
    )
    df.to_csv(f"{output}.csv", index_label="TAXID")

    with open(f"{output}.json", "w") as write_file:
        json.dump(res, write_file)
    logging.info(
        f"Step {7 if nb_steps == 7 else 8 }/{nb_steps}: writing sam2lca results:\n* JSON to {output}.json\n* CSV to {output}.csv"
    )

    return res


def get_db_connection(db_path):
    return rocksdb.DB(db_path, opts=rocksdb.Options(), read_only=True)
