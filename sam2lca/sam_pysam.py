#!/usr/bin/env python3

import pysam
from sam2lca.check_conserved_regions import (
    compute_coverage,
    flag_conserved_regions,
    is_in_conserved,
)
import rocksdb
from sam2lca.config import OPTS_read
from tqdm.contrib.concurrent import process_map
import logging
from tqdm import tqdm
from collections import ChainMap


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
        logging.error(e)
        return None


class Alignment:
    def __init__(self, al_file, nb_steps):
        """
        Args:
            al_file (str): Path to alignment file
            nb_steps (int): Number of steps in sam2lca
        """
        mode = {"sam": "r", "bam": "rb", "cram": "rc"}
        filetype = al_file.split(".")[-1]
        self.nb_steps = nb_steps
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

        logging.info(f"Step 2/{self.nb_steps}: Loading acc2tax database")
        DB = rocksdb.DB(dbname, opts=OPTS_read, read_only=True)

        logging.info(f"Step 3/{self.nb_steps}: Converting accession numbers to TAXIDs")
        self.acc2tax = accession_to_taxid_lookup(self.refs)

        del DB

        return self.acc2tax

    def get_conserved_regions(self, ref):
        """Get regions with higher than expected coverage, with zscore method

        Args:
            ref (pysam reference): one of pysam.alignment.references
        Returns:
            dict: {reference: conserved_regions[(start,end)]}
        """

        al_file = pysam.AlignmentFile(self.al_file, self.mode)
        refcov = compute_coverage(al_file.count_coverage(ref))
        window_size = min(al_file.get_reference_length(ref), 500)
        conserved_regions = flag_conserved_regions(refcov, window_size=window_size)
        return {ref: conserved_regions}

    def __get_reads_refs__(self, identity, minlength, check_conserved, process):
        """Get reads passing identity threshold for each reference

        Args:
            identity (float): identity threshold
            minlength(int): Length threshold.
            check_conserved(bool): Check if read is mapped in conserved region
            process(int): Number of processes
        """
        total_reads = int(pysam.view("-c", f"-@ {process}", self.al_file).rstrip())
        al_file = pysam.AlignmentFile(self.al_file, self.mode, threads=process)
        for read in tqdm(al_file, unit="reads", total=total_reads):
            if read.has_tag("NM") and not read.is_unmapped:
                mismatch = read.get_tag("NM")
                alnLen = read.query_alignment_length
                readLen = read.query_length
                ident = (alnLen - mismatch) / readLen
                if ident >= identity and alnLen >= minlength:
                    if check_conserved:
                        is_conserved = is_in_conserved(
                            read, self.cons_dict[read.reference_name]
                        )
                        if is_conserved is False:
                            try:
                                self.read_ref_dict.setdefault(
                                    read.query_name,
                                    set({self.acc2tax[read.reference_name]}),
                                ).add(self.acc2tax[read.reference_name])
                            except KeyError:
                                self.read_ref_dict.setdefault(
                                    read.query_name, set({0})
                                ).add(0)

                    else:
                        try:
                            self.read_ref_dict.setdefault(
                                read.query_name,
                                set({self.acc2tax[read.reference_name]}),
                            ).add(self.acc2tax[read.reference_name])
                        except KeyError:
                            self.read_ref_dict.setdefault(
                                read.query_name, set({0})
                            ).add(0)
        al_file.close()

    def get_reads(self, process=2, identity=0.9, minlength=30, check_conserved=False):
        """Get reads passing identity threshold
        Args:
            dbname(str): Path of RocksDB acc2tax database
            process (int, optional): Number of processes. Defaults to 2.
            identity (float, optional): Identity thresold. Defaults to 0.9.
            minlength(int, optional): Length threshold. Default to 30
            check_conserved(bool, optional): Check conserved regions
        Returns:
            dict: {read: set(mappeds TAXIDs)}
        """

        self.read_ref_dict = dict()

        if check_conserved:
            logging.info(f"Step 4/{self.nb_steps}: Getting conserved regions")

        if process == 1:
            self.cons_dict = dict()
            for ref in tqdm(self.refs):
                if check_conserved:
                    self.cons_dict.update(ref, self.get_conserved_regions(ref))

                self.__get_reads_single__(ref, identity, minlength, check_conserved)

        else:
            if check_conserved:
                cons_res = process_map(
                    self.get_conserved_regions,
                    self.refs,
                    max_workers=process,
                    chunksize=1,
                )
                self.cons_dict = dict(ChainMap(*cons_res))

            logging.info(
                f"Step {4 if self.nb_steps == 7 else 5 }/{self.nb_steps}: Parsing reads in alignment file"
            )

            self.__get_reads_refs__(
                identity=identity,
                minlength=minlength,
                check_conserved=check_conserved,
                process=process,
            )

        return self.read_ref_dict
