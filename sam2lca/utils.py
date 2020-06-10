import os


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
