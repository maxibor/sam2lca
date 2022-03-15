#!/usr/bin/env python3

import multiprocessing
from functools import partial
from sam2lca.config import NCBI
from ete3 import Tree
from tqdm.contrib.concurrent import process_map
import logging


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
            # print(list(map(str, taxids[1])))
            ancestor = tree.get_common_ancestor(list(map(str, tuple(taxids)))).name

        except ValueError as e:
            print(e)
    return {taxids: int(ancestor)}


def get_unique_taxids_combs(taxids_comb, lcas_dict, tree):
    """
    Args:
        taxids_comb: frozenset(taxids_combs)
        lcas_dict (dict): dict with unique combination of taxids as keys and taxids of LCA as values
    """
    lcas_dict.update(taxids_to_lca(taxids_comb, tree=tree))


def read_accession_to_taxid(read_name, read_dict, lcas_dict, ancestors_dict):
    """
    Args:
        read_name (str): read name
        read_dict (dict): dict with read names as keys and taxids as values
        lcas_dict (dict): dict with taxids as keys and taxids of LCA as values
        ancestors_dict (dict): dict with read names as keys and taxids of LCA as values
    """
    ancestors_dict.update({read_name: lcas_dict[frozenset(read_dict[read_name])]})


def compute_lca_multi(read_dict, tree, process):

    if tree:
        thetree = Tree(tree, format=1)
    else:
        thetree = None

    manager = multiprocessing.Manager()
    unique_taxids_combs = set({frozenset(i) for i in read_dict.values()})

    logging.info("Step 4/6: Getting unique TAXID combinations")
    if process == 1:
        lcas_dict = dict()  # dict(frozenset(taxids) : lca)
        ancestors_dict = dict()  # dict(read : lca)
        for taxids_combs in unique_taxids_combs:
            lcas_dict.update(taxids_to_lca(taxids_combs, thetree))
        logging.info("Step 5/6: Computing LCAs")
        for read in read_dict:
            print(lcas_dict[frozenset(read_dict[read])])
            ancestors_dict.update({read: lcas_dict[frozenset(read_dict[read])]})

    else:
        lcas_dict = manager.dict()  # dict(frozenset(taxids) : lca)
        ancestors_dict = manager.dict()

        get_unique_taxids_combs_partial = partial(
            get_unique_taxids_combs, lcas_dict=lcas_dict, tree=tree
        )
        process_map(
            get_unique_taxids_combs_partial,
            unique_taxids_combs,
            chunksize=1,
            max_workers=process,
        )

        logging.info("Step 5/6: Computing LCAs")
        read_accession_to_taxid_partial = partial(
            read_accession_to_taxid,
            read_dict=read_dict,
            lcas_dict=lcas_dict,
            ancestors_dict=ancestors_dict,
        )

        process_map(
            read_accession_to_taxid_partial,
            read_dict.keys(),
            chunksize=1,
            max_workers=process,
        )
    return ancestors_dict
