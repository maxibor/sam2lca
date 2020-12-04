#!/usr/bin/env python3

import multiprocessing
from functools import partial
import rocksdb
import ete3


class ReadToLca():
    def __init__(self, read, read_dict):
        """
        Args:
            read(str): read name
            read_dict (dict): Read as key, references accessions in list
        """
        self.read = read
        self.mapping = read_dict[read]

    def ref_to_taxid_single(self):
        """Get TAXID for a reference accession
        """
        # DB: global variable
        try:
            taxo_hits = [int(DB.get(bytes(i, encoding='utf8')))
                         for i in self.mapping]
            self.taxo_hits = list(set(taxo_hits))
        except TypeError as e:
            print(e)
            pass

    def compute(self, tree):
        if len(self.taxo_hits) == 1:
            ancestor = self.taxo_hits[0]
        else:
            if not tree:
                tree = NCBI.get_topology(
                    self.taxo_hits, intermediate_nodes=True)
            ancestor = (tree
                        .get_common_ancestor([str(i) for i in self.taxo_hits])
                        .name)
        return({self.read: int(ancestor)})


def compute_lca_read(read, read_dict, tree=None):
    r = ReadToLca(read, read_dict)
    r.ref_to_taxid_single()
    res = r.compute(tree)
    return(res)


def compute_lca_multi(read_dict, dbname, tree, update, process):
    global DB
    global NCBI
    NCBI = ete3.NCBITaxa()
    print("Loading Taxonomy database")
    DB = rocksdb.DB(dbname, opts=rocksdb.Options(), read_only=True)
    print("Finished loading Taxonomy database")

    if tree:
        thetree = ete3.Tree(tree, format=1)
    else:
        thetree = None

    compute_lca_partial = partial(
        compute_lca_read, read_dict=read_dict, tree=thetree)

    with multiprocessing.Pool(process) as p:
        allres = p.map(compute_lca_partial, list(read_dict.keys()))

    res = {}
    for r in allres:
        res.update(r)
    return(res)


if __name__ == "__main__":

    read_dict = {'shigella_dysenteriae_1483': ['NC_000913.3', 'NC_007607.1'],
                 'shigella_dysenteriae_1378': ['NC_007606.1', 'NC_007607.1']}
