import logging
from pysam import AlignmentFile, AlignedSegment
from tqdm import tqdm

def write_bam_by_taxid(
    target_taxids,
    infile,
    outfile_base,
    total_reads,
    read_taxid_dict,
    acc2tax_dict,
    taxid_info_dict,
    identity,
    edit_distance,
    minlength,
    process
):
    """Write alignment file with taxonomic information tag

    Args:
        target_taxids (list): list of taxid of target taxons
        infile (str): Path to alignment file
        outfile (str): Path to output file
        total_reads(int): Total number of reads in file
        read_taxid_dict (dict): Dictionary of read with their corresponding LCA taxid
        taxid_info_dict(dict): Dictionary of taxid with their count and metadata
        acc2tax_dict(dict): Accession to taxid dictionary
        identity (float): Minimum identity for alignment
        edit_distance (int): Minimum edit distance for alignment
        minlength (int): Minimum length for alignment
        process (int): Number of threads to use
    """
    mode = {"sam": "r", "bam": "rb", "cram": "rc"}
    filetype = infile.split(".")[-1]

    with AlignmentFile(infile, mode[filetype], threads=process) as samfile:
        aligned_refs = set()
        for ref_stat in samfile.get_index_statistics():
            refname = ref_stat[0]
            nb_mapped_reads = ref_stat[1]
            if nb_mapped_reads > 0:
                aligned_refs.add(refname)
        
        reflens_dict = dict()
        for taxid in target_taxids:
            reflens_dict[taxid] = dict()
            for ref in aligned_refs:
                try:
                    if taxid in acc2tax_dict[ref].taxid_lineage:
                        reflens_dict[taxid][ref] = samfile.get_reference_length(ref)
                except KeyError:
                    logging.error(f"Error with reference {ref}")
                    continue
        header_dict = dict()        
        for taxid in target_taxids:
            header_dict[taxid] = {
                "HD": {"VN": "1.0", "SO": "unsorted"},
                "SQ": [
                    {"LN": reflens_dict[taxid][refname], "SN": refname} for refname in reflens_dict[taxid].keys()
                ]
            }
        refid_dict = dict()
        for taxid in target_taxids:
            refid_dict[taxid] = dict()
            for i, ref in enumerate(reflens_dict[taxid].keys()):
                refid_dict[taxid][ref] = i 

        open_files_dict = dict()
        for taxid in target_taxids:
            outfile = f"{outfile_base.replace('.sam2lca.bam','')}_taxid_{taxid}.sam2lca.bam"
            logging.info(f"  - {taxid_info_dict[taxid]['name']} (taxid: {taxid}) - {outfile}")
            open_files_dict[taxid] = AlignmentFile(
                outfile, 
                "wb", 
                header=header_dict[taxid], 
                reference_names=refid_dict[taxid].keys(), threads=process)

        for read in tqdm(samfile, unit="reads", total=total_reads):
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
                    read_taxid = read_taxid_dict[read.query_name]
                    taxid_intersection = list(set(target_taxids).intersection(set(taxid_info_dict[read_taxid]['lineage_taxid'])))
                    if len(taxid_intersection) == 1 :
                        target_taxid = taxid_intersection[0]
                        try:
                            rw = AlignedSegment()
                            rw.query_name = read.query_name
                            rw.query_sequence= read.query_sequence
                            rw.flag = read.flag
                            rw.reference_id = refid_dict[target_taxid][read.reference_name]
                            rw.reference_start = read.reference_start
                            rw.mapping_quality = read.mapping_quality
                            rw.cigar = read.cigar
                            if read.next_reference_name:
                                rw.next_reference_id = refid_dict[target_taxid][read.next_reference_name]
                                rw.next_reference_start=read.next_reference_start
                            rw.template_length=read.template_length
                            rw.query_qualities = read.query_qualities
                            rw.tags = read.tags

                            rw.set_tag("XT", read_taxid, "i")
                            rw.set_tag(
                                "XN",
                                taxid_info_dict[read_taxid]["name"],
                                "Z",
                            )
                            rw.set_tag(
                                "XR",
                                taxid_info_dict[read_taxid]["rank"],
                                "Z",
                            )
                            open_files_dict[target_taxid].write(rw)
                        except KeyError as e:
                            logging.error(e)
                            continue
        for taxid in target_taxids:
            open_files_dict[taxid].close()

                                


def write_bam_tags(
    infile,
    outfile,
    total_reads,
    read_taxid_dict,
    taxid_info_dict,
    identity,
    edit_distance,
    minlength,
    process,
    unclassified_taxid=12908,
    
):
    """Write alignment file with taxonomic information tag

    Args:
        infile (str): Path to alignment file
        outfile (str): Path to output file
        total_reads(int): Total number of reads in file
        read_taxid_dict (dict): Dictionary of read with their corresponding LCA taxid
        taxid_info_dict(dict): Dictionary of taxid with their count and metadata
        identity (float): Minimum identity for alignment
        edit_distance (int): Minimum edit distance for alignment
        minlength (int): Minimum length for alignment
        process (int): Number of threads to use
    """
    logging.info(f"* Writing BAM file with taxonomic information tag to {outfile}")
    mode = {"sam": "r", "bam": "rb", "cram": "rc"}
    filetype = infile.split(".")[-1]
    with AlignmentFile(infile, mode[filetype], threads=process) as samfile:
        with AlignmentFile(outfile, "wb", template=samfile, threads=process) as bamout:
            for read in tqdm(samfile, unit="reads", total=total_reads):
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
                        try:
                            taxid = read_taxid_dict[read.query_name]
                            read.set_tag("XT", taxid, "i")
                            read.set_tag(
                                "XN",
                                taxid_info_dict[taxid]["name"],
                                "Z",
                            )
                            read.set_tag(
                                "XR",
                                taxid_info_dict[taxid]["rank"],
                                "Z",
                            )
                        except KeyError:
                            read.set_tag("XT", unclassified_taxid, "i")
                            read.set_tag(
                                "XN",
                                "unclassified sequences",
                                "Z",
                            )
                    else:
                        read.set_tag("XT", unclassified_taxid, "i")
                        read.set_tag(
                            "XN",
                            "unclassified sequences",
                            "Z",
                        )
                else:
                    read.set_tag("XT", unclassified_taxid, "i")
                    read.set_tag(
                        "XN",
                        "unclassified sequences",
                        "Z",
                    )
                bamout.write(read)
