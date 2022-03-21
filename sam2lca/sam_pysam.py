#!/usr/bin/env python3

import pysam
from functools import partial
from multiprocessing import Manager
from sam2lca.check_conserved_regions import (
    compute_coverage,
    flag_conserved_regions,
    is_in_conserved,
)
import rocksdb
from sam2lca.config import OPTS_read
from tqdm.contrib.concurrent import thread_map, process_map
import logging
from tqdm import tqdm


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
        present_refs = set()
        for ref_stat in alignment.get_index_statistics():
            refname = ref_stat[0]
            nb_mapped_reads = ref_stat[1]
            if nb_mapped_reads > 0:
                present_refs.add(refname)
        self.refs = tuple(present_refs)
        alignment.close()

    def get_refs_taxid(self, dbname):
        """Get taxids of referernce sequences

        Args:
            dbname (str): Path of RocksDB acc2tax database
        Returns:
            dict: {ref: taxid}
        """

        global DB

        logging.info("Step 1/6: Loading acc2tax database")
        DB = rocksdb.DB(dbname, opts=OPTS_read, read_only=True)

        logging.info("Step 2/6: Converting accession numbers to TAXIDs")
        self.acc2tax = accession_to_taxid_lookup(self.refs)

        del DB

        return self.acc2tax

    def __get_reads_single__(
        self, ref, identity, minlength, check_conserved, read_ref_dict
    ):
        """Get reads passing identity threshold for each reference

        Args:
            ref (pysam reference): one of pysam.alignment.references
            identity (float): identity threshold
            minlength(int): Length threshold.
            check_conserved(bool): Check if read is mapped in conserved region
            read_ref_dict(dict): {read:set(taxid of mapped_references)}
        """
        al_file = pysam.AlignmentFile(self.al_file, self.mode)
        refcov = compute_coverage(al_file.count_coverage(ref))
        window_size = min(al_file.get_reference_length(ref), 500)
        if check_conserved:
            conserved_regions = flag_conserved_regions(refcov, window_size=window_size)
        reads = al_file.fetch(ref)
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
                            try:
                                read_ref_dict.setdefault(
                                    read.query_name,
                                    set({self.acc2tax[read.reference_name]}),
                                ).add(self.acc2tax[read.reference_name])
                            except KeyError:
                                read_ref_dict.setdefault(read.query_name, set({0})).add(
                                    0
                                )

                    else:
                        try:
                            read_ref_dict.setdefault(
                                read.query_name,
                                set({self.acc2tax[read.reference_name]}),
                            ).add(self.acc2tax[read.reference_name])
                        except KeyError:
                            read_ref_dict.setdefault(read.query_name, set({0})).add(0)
        al_file.close()

    def get_reads(
        self,
        process=2,
        identity=0.9,
        minlength=30,
        check_conserved=False,
    ):
        """Get reads passing identity threshold
        Args:
            dbname(str): Path of RocksDB acc2tax database
            process (int, optional): Number of processes. Defaults to 2.
            identity (float, optional): Identity thresold. Defaults to 0.9.
            minlength(int, optional): Length threshold. Default to 30
            check_conserved(bool, optional): Check conserved regions

        """

        read_ref_dict = dict()
        if process == 1:
            for ref in self.refs:
                self.__get_reads_single__(
                    ref,
                    identity,
                    minlength,
                    check_conserved,
                    read_ref_dict,
                )
        else:
            get_reads_partial = partial(
                self.__get_reads_single__,
                identity=identity,
                minlength=minlength,
                check_conserved=check_conserved,
                read_ref_dict=read_ref_dict,
            )
            logging.info("Step 3/6: Parsing reads in alignment file")
            thread_map(get_reads_partial, self.refs, max_workers=process, chunksize=1)

        return read_ref_dict
