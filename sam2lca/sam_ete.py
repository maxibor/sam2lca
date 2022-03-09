#!/usr/bin/env python3

import multiprocessing
from functools import partial
from sam2lca.config import NCBI, OPTS_read
from ete3 import Tree
import rocksdb
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from ordered_set import OrderedSet
import logging


def accession_to_taxid_lookup(accession_list):
    """Get TAXID of hits from hit accessions per read

    Args:
        accession list (list): list of unique accession numbers
    Returns:
        dict: {accession: taxid}
    """
    try:
        return {i: int(DB.get(bytes(i, encoding="utf8"))) for i in tqdm(accession_list)}
    except TypeError as e:
        print(e)
        return None


def read_accession_to_taxid(read_hit_record, accession_taxid_dict):
    """Get TAXID of hits from hit accessions per read

    Args:
        read_hit_record (tuple): (read_name(str), accession(list))
    Returns:
        dict: {read_name(str): (tuple of taxids)}
    """
    try:
        return {
            read_hit_record[0]: tuple(
                OrderedSet([accession_taxid_dict[i] for i in read_hit_record[1]])
            )
        }
    except TypeError as e:
        print(e)
        return None


def taxids_to_lca(taxids, tree):
    """Run LCA on list of TAXID

    Args:
        taxids (set): tuple of taxids
        tree (ete tree): taxid tree
    """
    taxids = taxids
    if len(taxids) == 1:
        ancestor = taxids[0]
    else:
        if not tree:
            try:
                tree = NCBI.get_topology(taxids, intermediate_nodes=True)
            except ValueError as e:
                print(e)
        try:
            # print(list(map(str, taxids[1])))
            ancestor = tree.get_common_ancestor(list(map(str, taxids))).name

        except ValueError as e:
            print(e)
    return {taxids: int(ancestor)}


def compute_lca_multi(read_dict, accession_list, dbname, tree, process):
    global DB

    logging.info("Step 2/6: Loading Taxonomy database")

    DB = rocksdb.DB(dbname, opts=OPTS_read, read_only=True)
    logging.info("* Finished loading Taxonomy database")
    logging.info("Step 3/6: Converting accession numbers to TAXIDs")
    accession2taxid = accession_to_taxid_lookup(accession_list)
    del DB

    if tree:
        thetree = Tree(tree, format=1)
    else:
        thetree = None

    # if process >= 2:
    #     compute_lca_partial = partial(
    #         compute_lca_read, read_dict=read_dict, tree=thetree)

    #     with multiprocessing.Pool(process) as p:
    #         allres = p.map(compute_lca_partial, list(read_dict.keys()))
    # res = {}
    # for r in allres:
    #     res.update(r)
    # return(res)
    # else:

    ancestors = {}
    logging.info("Step 4/6: Getting unique TAXID combinations")
    if process == 1:
        for read in tqdm(read_dict.items()):
            # Getting TAXIDs the of accesions mapped to each read
            taxid_hits = read_accession_to_taxid(read, accession2taxid)
            logging.info("Step 5/6: Computing LCA")
            if taxid_hits:
                # Getting the LCA of these TAXIDs
                ancestors.update(taxids_to_lca(taxid_hits[read[0]], thetree))

    else:
        taxid_hits = {}
        read_accession_to_taxid_partial = partial(
            read_accession_to_taxid, accession_taxid_dict=accession2taxid
        )
        taxid_res = process_map(
            read_accession_to_taxid_partial,
            read_dict.items(),
            chunksize=1,
            max_workers=process,
        )
        # Creating dictornary with read name as key and tuple of TAXIDs as value
        [taxid_hits.update(r) for r in taxid_res if r]

        # Summarise reads with identical taxid combinations
        taxid_combs = (
            {}
        )  # taxids as keys, reads having mapping to these taxids as values
        for read, taxids in taxid_hits.items():
            if taxids not in taxid_combs:
                taxid_combs[taxids] = [read]
            else:
                taxid_combs[taxids].append(read)

        # Infer ancestors on unique combinations
        logging.info(f"* {len(taxid_combs)} unique combinations of taxa were found.")
        lca = {}
        logging.info("Step 5/6: Computing LCA")
        for taxid_res in tqdm(taxid_combs):  # taxids as key, taxid of LCA as value
            lca.update(taxids_to_lca(taxid_res, thetree))
        # Extract for each read ancestor taxid
        for taxids, reads in taxid_combs.items():
            for read in reads:
                ancestors.update({read: lca[taxids]})

    return ancestors


if __name__ == "__main__":

    read_dict = {
        "shigella_dysenteriae_1483": ["NC_000913.3", "NC_007607.1"],
        "shigella_dysenteriae_1378": ["NC_007606.1", "NC_007607.1"],
    }
