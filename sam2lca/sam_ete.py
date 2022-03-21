#!/usr/bin/env python3

import multiprocessing
from functools import partial
from sam2lca.config import NCBI
from ete3 import Tree
from tqdm.contrib.concurrent import process_map, thread_map
import logging
import sys
from collections import ChainMap


def taxids_to_lca(taxids, tree):
    """Run LCA on list of TAXID

    Args:
        taxids (frozenset): frozenset of taxids
        tree (ete tree): taxid tree
    """

    if len(taxids) == 1:
        ancestor = tuple(taxids)[0]
    else:
        if not tree:
            try:
                tree = NCBI.get_topology(tuple(taxids), intermediate_nodes=True)
            except ValueError as e:
                print(e)
        try:
            ancestor = tree.get_common_ancestor(list(map(str, tuple(taxids)))).name

        except ValueError as e:
            print(e)
    return {taxids: int(ancestor)}


def read_accession_to_taxid(read_name, read_dict, lcas_dict):
    """
    Args:
        read_name (str): read name
        read_dict (dict): dict with read names as keys and taxids as values
        lcas_dict (dict): dict with taxids as keys and taxids of LCA as values
    """
    ancestors_dict.update({read_name: lcas_dict[frozenset(read_dict[read_name])]})


def compute_lca_multi(read_dict, tree, process):

    if tree:
        thetree = Tree(tree, format=1)
    else:
        thetree = None

    global ancestors_dict
    ancestors_dict = dict()  # dict(read : lca)

    unique_taxids_combs = set({frozenset(i) for i in read_dict.values()})
    i = 0
    for t in unique_taxids_combs:
        i += 1
        if i == 10:
            break

    logging.info("Step 4/6: Getting unique TAXID combinations")
    lcas_dict = dict()
    if process == 1:
        for taxids_combs in unique_taxids_combs:
            lcas_dict.update(taxids_to_lca(taxids_combs, thetree))
        logging.info("Step 5/6: Computing LCAs")
        for read in read_dict:
            ancestors_dict.update({read: lcas_dict[frozenset(read_dict[read])]})

    else:
        get_taxids_partial = partial(taxids_to_lca, tree=tree)
        lcas_res = process_map(
            get_taxids_partial,
            unique_taxids_combs,
            chunksize=1,
            max_workers=process,
        )

        lcas_dict = dict(ChainMap(*lcas_res))

        logging.info("Step 5/6: Assigning LCA to reads")
        read_accession_to_taxid_partial = partial(
            read_accession_to_taxid, read_dict=read_dict, lcas_dict=lcas_dict
        )

        thread_map(
            read_accession_to_taxid_partial,
            read_dict.keys(),
            chunksize=1,
            max_workers=process,
        )

    return ancestors_dict
