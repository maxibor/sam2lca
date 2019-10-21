#!/usr/bin/env python3

import pysam
from functools import partial
import multiprocessing

class alignment():
    def __init__(self, al_file, filetype):
        """
        Args:
            al_file (str): Path to alignment file
            filetype (str): Type of alignment file (sam, bam, cram)
        """
        mode = {'sam':'r', 'bam':'rb', 'cram':'rc'}
        alignment = pysam.AlignmentFile(al_file, mode[filetype])
        self.al_file = al_file
        self.mode = mode[filetype]
        self.refs = self.alignment.references
        alignment.close()

    def __get_reads_single__(self, ref, identity):
        """Get reads passing identity threshold for each reference
        
        Args:
            ref (pysam reference): one of pysam.alignment.references
            identity (float): identity threshold
        """
        resdic = {}
        al_file = pysam.AlignmentFile(self.al_file, self.mode)
        reads = al_file.fetch(ref, multiple_iterators=True)
        for read in reads:
            mismatch = read.get_tag("NM")
            alnLen = read.query_alignment_length
            readLen = read.query_length
            ident = (alnLen - mismatch) / readLen
            if ident >= identity:
                resdic[read.query_name] = reads.reference_name
        return(resdic)


    def get_reads(self, process=2, identity=0.9):
        """Get reads passing identity threshold
        
        Args:
            process (int, optional): Number of processes. Defaults to 2.
            identity (float, optional): Identity thresold. Defaults to 0.9.
        """
        get_reads_partial = partial(__get_reads_single__, identity=identity)
        with multiprocessing.Pool(process) as p:
            results = p.map(get_reads_partial, self.refs)
        
        def merge_dict(dict_list):
            """Merge list of dictionaries by key while preserving
            values
            
            Args:
                dict_list (list of dictionaries)
            """
            resdict = {}
            alldict = [list(i.items())[0] for i in dict_list]
            for i in alldict:
                if i[0] not in alldict.keys():
                    resdict[i[0]] = [i[1]]
                else:
                    resdict[i[0]].append(i[1])
            return(resdict)
        
        self.mapping = merge_dict(results)