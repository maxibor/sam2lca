#!/usr/bin/env python3

import multiprocessing
from functools import partial
from sam2lca.config import NCBI, OPTS_read
from ete3 import Tree
import rocksdb
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

def accession_hits_to_taxid(read_hit_record):
    """Get TAXID of hits from hit accessions per read

    Args:
        read_hit_record (tuple): (read_name(str), accession(list))
    """
    try:
        return {read_hit_record[0]: [int(DB.get(bytes(i, encoding='utf8')))
                        for i in read_hit_record[1]]}
    except TypeError as e:
        print(e)
        return None   

def taxids_to_lca(taxids, tree):
    """Run LCA on list of TAXID

    Args:
        taxids (list): list of taxids
        tree (ete tree): taxid tree
    """
    if len(list(set(taxids))) == 1:
        ancestor = taxids[0]
    else:
        if not tree:
            try:
                tree = NCBI.get_topology(
                    taxids, intermediate_nodes=True)
            except ValueError as e:
                print(e)
                print(taxids)
        try:
            ancestor = (tree
                        .get_common_ancestor([str(i) for i in taxids])
                        .name)
        except ValueError as e:
            print(e)
            print(taxids)
    return({",".join([str(t) for t in taxids]): int(ancestor)})


def compute_lca_multi(read_dict, dbname, tree, process):
    global DB

    print("Loading Taxonomy database")

    DB = rocksdb.DB(dbname, opts=OPTS_read, read_only=True)
    print("Finished loading Taxonomy database")

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
    print("Computing LCA")
    if process == 1:
        for read in tqdm(read_dict.items()):
            taxid_hits = accession_hits_to_taxid(read)

            del(DB)
            
            if taxid_hits:
                ancestors.update(taxids_to_lca((read[0], taxid_hits[read[0]]), thetree))
            
    else:
        taxid_hits = {}
        taxid_res = process_map(accession_hits_to_taxid, read_dict.items(), chunksize=1, max_workers=process)
        [taxid_hits.update(r) for r in taxid_res if r]

        del(DB)

        # Summarise reads with identical taxid combinations
        taxid_combs = {}
        for read, taxids in taxid_hits.items():
            taxids_hash = ",".join([str(t) for t in taxids])
            if taxids_hash not in taxid_combs:
                taxid_combs[taxids_hash] = [read]
            else:
                taxid_combs[taxids_hash].append(read)

        # Infer ancestors on unique combinations
        unique_taxid_hashs = list(set(taxid_combs.keys()))
        print(f"{len(unique_taxid_hashs)} unique combinations of taxa were found.")
        lca = {}
        for taxid_res in tqdm(unique_taxid_hashs):
            lca.update(taxids_to_lca([int(t) for t in taxid_res.split(",")], thetree))

        # Extract for each read ancestor taxid
        for taxids_hash, reads in taxid_combs.items():
            for read in reads:
                ancestors.update({read: lca[taxids_hash]})

    return(ancestors)


if __name__ == "__main__":

    read_dict = {'shigella_dysenteriae_1483': ['NC_000913.3', 'NC_007607.1'],
                 'shigella_dysenteriae_1378': ['NC_007606.1', 'NC_007607.1']}
