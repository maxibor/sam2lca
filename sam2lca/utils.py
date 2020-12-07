import os
import json
import rocksdb
from sam2lca.config import NCBI
import pandas as pd

def get_script_dir():
    return(os.path.dirname(os.path.realpath(__file__)))


def count_reads_taxid(read_taxid_dict):
    """Returns number of reads matching TAXID

    Args:
        read_taxid_dict (dict): {read_name(str): TAXID(int)}
    """
    taxid_cnt = {}
    for r in read_taxid_dict:
        if read_taxid_dict[r] not in taxid_cnt:
            taxid_cnt[read_taxid_dict[r]] = 1
            continue
        else:
            taxid_cnt[read_taxid_dict[r]] += 1
    return(taxid_cnt)


def output_file(sam_path):
    out = os.path.basename(sam_path).split(".")[:-1]
    out = ".".join(out)+".sam2lca"
    return(out)


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
    modes = {'bam': 'rb', 'sam': 'r', 'cram': 'rc'}
    try:
        return(modes[extension])
    except KeyError:
        raise Exception(f"{extension} file extension not supported")


def taxid_to_lineage(taxid_count_dict, output):
    res = {}
    for taxid in taxid_count_dict:
        read_count = taxid_count_dict[taxid]
        sciname = NCBI.get_taxid_translator([taxid])[taxid]
        rank = NCBI.get_rank([taxid])[taxid]
        taxid_lineage = NCBI.get_lineage(taxid)
        scinames = NCBI.get_taxid_translator(taxid_lineage)
        ranks = NCBI.get_rank(taxid_lineage)
        lineage = [{ranks[taxid]:scinames[taxid]} for taxid in taxid_lineage]
        res[taxid] = {'name': sciname,
                      'rank': rank,
                      'count': read_count,
                      'lineage': lineage}

    (pd.DataFrame(res).transpose().sort_values("count", ascending=False)
    .to_csv(f"{output}.csv", index_label='TAXID'))
    with open(f"{output}.json", "w") as write_file:
        json.dump(res, write_file)
    print(f"sam2lca results written to:\n- {output}.json\n- {output}.csv")

    return(res)

def get_db_connection(db_path):
    return(rocksdb.DB(db_path, opts=rocksdb.Options(), read_only=True))
