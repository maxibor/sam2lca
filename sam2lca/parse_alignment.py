#!/usr/bin/env python3

import pysam
from sam2lca.check_conserved_regions import (
    compute_coverage,
    flag_conserved_regions,
    is_in_conserved,
)
import rocksdb
from sam2lca.rocksdb import OPTS_read
from tqdm.contrib.concurrent import process_map
import logging
from tqdm import tqdm
from collections import ChainMap
from taxopy import Taxon
from taxopy.exceptions import TaxidError


def accession_to_taxid_lookup(accession_list, taxo_db, unclassified_taxid=12908):
    """Get TAXID of hits from hit accessions per read

    Args:
        accession list (list): list of unique accession numbers
        taxo_db (taxopy database): taxonomy database
        unclassified_taxid(int): TAXID for unclassified sequences
    Returns:
        dict: {accession: taxid}
    """
    acc2tax = dict()
    for i in tqdm(accession_list):
        try:
            taxid = int(DB.get(bytes(i, encoding="utf8")))
            if taxid != 0:
                acc2tax[i] = Taxon(taxid, taxo_db)
        except (TypeError, TaxidError) as e:
            acc2tax[i] = Taxon(unclassified_taxid, taxo_db)
            continue
    return acc2tax


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

    def get_refs_taxid(self, acc2tax_dbname, taxo_db):
        """Get taxids of reference sequences

        Args:
            acc2tax_dbname (str): Path of RocksDB acc2tax database
            taxo_db (taxopy database): taxonomy database
        Returns:
            dict: {ref: taxid}
        """

        global DB

        logging.info(f"Step 2/{self.nb_steps}: Loading acc2tax database")
        DB = rocksdb.DB(acc2tax_dbname, opts=OPTS_read, read_only=True)

        logging.info(f"Step 3/{self.nb_steps}: Converting accession numbers to TAXIDs")
        self.acc2tax = accession_to_taxid_lookup(self.refs, taxo_db)

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

    def __get_reads_refs__(
        self, identity, edit_distance, minlength, check_conserved, process, unclassified_taxid=12908
    ):
        """Get reads passing identity threshold

        Args:
            identity (float): identity threshold
            edit_distance (int): edit distance threshold
            minlength(int): Length threshold.
            check_conserved(bool): Check if read is mapped in conserved region
            process(int): Number of processes
        """
        self.total_reads = int(pysam.view("-c", f"-@ {process}", self.al_file).rstrip())
        al_file = pysam.AlignmentFile(self.al_file, self.mode, threads=process)
        for read in tqdm(al_file, unit="reads", total=self.total_reads):
            if read.has_tag("NM") and not read.is_unmapped:
                mismatch = read.get_tag("NM")
                alnLen = read.query_alignment_length
                readLen = read.query_length
                if edit_distance:
                    threshold = edit_distance
                    align_value = mismatch
                else:
                    threshold = 1 - identity
                    align_value = 1 - ((alnLen - mismatch) / readLen)
                if align_value <= threshold and alnLen >= minlength:
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
                                    read.query_name, set({unclassified_taxid})
                                ).add(unclassified_taxid)

                    else:
                        try:
                            self.read_ref_dict.setdefault(
                                read.query_name,
                                set({self.acc2tax[read.reference_name]}),
                            ).add(self.acc2tax[read.reference_name])
                        except KeyError:
                            self.read_ref_dict.setdefault(
                                read.query_name, set({unclassified_taxid})
                            ).add(unclassified_taxid)
        al_file.close()

    def get_reads(self, process=2, identity=0.9, edit_distance=None, minlength=30, check_conserved=False):
        """Get reads passing identity threshold
        Args:
            dbname(str): Path of RocksDB acc2tax database
            process (int, optional): Number of processes. Defaults to 2.
            identity (float, optional): Identity thresold. Defaults to 0.9.
            edit_distance (int): edit distance threshold
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
                    self.cons_dict.update(ref, self.get_conserved_regions(ref))
            else:
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
            edit_distance=edit_distance,
            minlength=minlength,
            check_conserved=check_conserved,
            process=process,
        )

        return self.read_ref_dict
