import os
import json
import rocksdb
from taxopy import Taxon
import pandas as pd
import logging
import warnings
from collections import Counter, ChainMap
from tqdm.contrib.concurrent import thread_map
from functools import partial

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
    except Exception:
        sciname = "unclassified sequences"
        rank = "NA"
        lineage = "NA"

    taxid_info_dict[taxid] = {
        "name": sciname,
        "rank": rank,
        "count_taxon": read_count_taxon,
        "count_descendant": read_count_descendant,
        "lineage": lineage,
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


def get_db_connection(db_path):
    return rocksdb.DB(db_path, opts=rocksdb.Options(), read_only=True)
