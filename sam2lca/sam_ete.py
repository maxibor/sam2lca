#!/usr/bin/env python3

import multiprocessing
from functools import partial
from sam2lca.config import NCBI, OPTS_read
from ete3 import Tree
import rocksdb
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool

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

def taxids_to_lca(read_record, tree):
    """Run LCA on list of TAXID

    Args:
        read_record(tuple): (read_name(str), taxid_hist(list))
        tree (ete tree): taxid tree
    """    
    read = read_record[0]
    taxids = read_record[1]
    if len(taxids) == 1:
        ancestor = taxids[0]
    else:
        if not tree:
            tree = NCBI.get_topology(
                taxids, intermediate_nodes=True)
        ancestor = (tree
                    .get_common_ancestor([str(i) for i in taxids])
                    .name)
    return({read:int(ancestor)})


def compute_lca_multi(read_dict, dbname, tree, update, process):
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
        with ThreadPool(process) as p:
            taxid_res = p.map(accession_hits_to_taxid, read_dict.items())
        [taxid_hits.update(r) for r in taxid_res if r]
        
        del(DB)

        taxids_to_lca_multi = partial(taxids_to_lca, tree=thetree)
        with multiprocessing.Pool(process) as p:
            res = p.map(taxids_to_lca_multi, taxid_hits.items())
        [ancestors.update(r) for r in res]

    return(ancestors)
    

   


if __name__ == "__main__":

    read_dict = {'shigella_dysenteriae_1483': ['NC_000913.3', 'NC_007607.1'],
                 'shigella_dysenteriae_1378': ['NC_007606.1', 'NC_007607.1']}
