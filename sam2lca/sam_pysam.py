#!/usr/bin/env python3

import pysam
from functools import partial
from sam2lca.check_conserved_regions import (
    compute_coverage,
    flag_conserved_regions,
    is_in_conserved,
)
import multiprocessing


class Alignment:
    def __init__(self, al_file):
        """
        Args:
            al_file (str): Path to alignment file
        """
        mode = {"sam": "r", "bam": "rb", "cram": "rc"}
        filetype = al_file.split(".")[-1]
        alignment = pysam.AlignmentFile(al_file, mode[filetype])
        self.al_file = al_file
        self.mode = mode[filetype]
        self.refs = alignment.references
        alignment.close()

    def __get_reads_single__(self, ref, identity, minlength, check_conserved):
        """Get reads passing identity threshold for each reference

        Args:
            ref (pysam reference): one of pysam.alignment.references
            identity (float): identity threshold
            minlength(int): Length threshold.
            check_conserved(bool): Check if read is mapped in conserved region
        """
        resdic = {}
        al_file = pysam.AlignmentFile(self.al_file, self.mode)
        refcov = compute_coverage(al_file.count_coverage(ref))
        window_size = min(al_file.get_reference_length(ref), 500)
        if check_conserved:
            conserved_regions = flag_conserved_regions(refcov, window_size=window_size)
        reads = al_file.fetch(ref, multiple_iterators=True)
        for read in reads:
            if read.has_tag("NM"):
                mismatch = read.get_tag("NM")
                alnLen = read.query_alignment_length
                readLen = read.query_length
                ident = (alnLen - mismatch) / readLen
                if ident >= identity and alnLen >= minlength:
                    if check_conserved:
                        is_conserved = is_in_conserved(read, conserved_regions)
                        if is_conserved is False:
                            resdic[read.query_name] = read.reference_name
                    else:
                        resdic[read.query_name] = read.reference_name
        return resdic

    def get_reads(self, process=2, identity=0.9, minlength=30, check_conserved=False):
        """Get reads passing identity threshold

        Args:
            process (int, optional): Number of processes. Defaults to 2.
            identity (float, optional): Identity thresold. Defaults to 0.9.
            minlength(int, optional): Length threshold. Default to 30
            check_conserved(bool, optional): Check conserved regions
        """
        get_reads_partial = partial(
            self.__get_reads_single__,
            identity=identity,
            minlength=minlength,
            check_conserved=check_conserved,
        )

        #####################################
        ## Debugging non-parallelized loop ##
        #####################################
        results = []
        for ref in self.refs:
            results.append(
                self.__get_reads_single__(
                    ref=ref,
                    identity=identity,
                    minlength=minlength,
                    check_conserved=check_conserved,
                )
            )

        # with multiprocessing.Pool(process) as p:
        #     results = p.map(get_reads_partial, self.refs)

        def merge_dict(dict_list):
            """Merge list of dictionaries by key while preserving
            values

            Args:
                dict_list (list of dictionaries)
            Returns:
                dict: {read:[mapped_references]}
            """
            resdict = {}
            alldict = [list(i.items()) for i in dict_list]
            for r in alldict:
                for i in r:
                    if i[0] not in resdict.keys():
                        resdict[i[0]] = [i[1]]
                    else:
                        resdict[i[0]].append(i[1])
            return resdict

        return merge_dict(results)
