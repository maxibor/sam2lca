#!/usr/bin/env python3

import ete3
import multiprocessing
from functools import partial


class lca():
    def __init__(self, mapping, tree=None):
        """
        Args:
            mapping (dict of list): Reads as keys, references in list
            tree (str): Path to Newick phylogenetic/Taxonomic tree. Defaults to None.
        """
        self.mapping = mapping
        if tree:
            self.tree = ete3.Tree(tree)
        else:
            self.tree = None
    
    def compute(self, process):
        result = {}
        if not self.tree:
            ncbi = ete3.NCBITaxa()
            ncbi.update_taxonomy_database()

        def get_lca(key):
            resdict = {}
            if len(self.maping[key] > 1):
                ancestor = str(self.maping[key])
                
            else:
                if not self.tree:
                    tree = ncbi.get_topology(self.maping[key])
                else:
                    tree = self.tree
                ancestor = (tree
                            .get_common_ancestor
                            ([str(i) for i in self.maping[key]])
                            .name)

            return(resdict[key] = ancestor)
        
        with multiprocessing.Pool(process) as p:
            res = p.map(self.mapping.keys())

        for l in res:
            result = result.update(l)
        return(result)


