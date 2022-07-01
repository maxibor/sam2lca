import logging
from pysam import AlignmentFile
from tqdm import tqdm

def write_bam_by_taxid(
    target_taxid,
    infile,
    outfile,
    total_reads,
    read_taxid_dict,
    acc2tax_dict,
    taxid_info_dict,
    identity,
    edit_distance,
    minlength
):
    """Write alignment file with taxonomic information tag

    Args:
        target_taxid (int): Taxid of target taxon
        infile (str): Path to alignment file
        outfile (str): Path to output file
        total_reads(int): Total number of reads in file
        read_taxid_dict (dict): Dictionary of read with their corresponding LCA taxid
        taxid_info_dict(dict): Dictionary of taxid with their count and metadata
        acc2tax_dict(dict): Accession to taxid dictionary
        identity (float): Minimum identity for alignment
        edit_distance (int): Minimum edit distance for alignment
        minlength (int): Minimum length for alignment
    """
    mode = {"sam": "r", "bam": "rb", "cram": "rc"}
    filetype = infile.split(".")[-1]
    outfile = f"{outfile.replace('.sam2lca.bam','')}_taxid_{target_taxid}.sam2lca.bam"

    with AlignmentFile(infile, mode[filetype]) as samfile:
        reflens = dict()

        for ref in samfile.references:
            try:
                if target_taxid in acc2tax_dict[ref].taxid_lineage:
                    reflens[ref] = samfile.get_reference_length(ref)
            except KeyError:
                continue

        header = {
            "HD": {"VN": "1.0", "SO": "unsorted"},
            "SQ": [
                {"LN": reflens[refname], "SN": refname} for refname in reflens.keys()
            ]
        }
        refid_dict = dict()
        for i, ref in enumerate(reflens.keys()):
            refid_dict[ref] = i 
        with AlignmentFile(outfile, "wb", header=header) as outfile:
            for read in tqdm(samfile,  desc=str(target_taxid), unit="reads", total=total_reads):
                if read.has_tag("NM") and not read.is_unmapped:
                    mismatch = read.get_tag("NM")
                    alnLen = read.query_alignment_length
                    readLen = read.query_length
                    ident = (alnLen - mismatch) / readLen
                    if edit_distance:
                        threshold = edit_distance
                        align_value = mismatch
                    else:
                        threshold = 1 - identity
                        align_value = 1 - ((alnLen - mismatch) / readLen)
                    if align_value <= threshold and alnLen >= minlength:
                        try:
                            taxid = read_taxid_dict[read.query_name]
                            if target_taxid in taxid_info_dict[taxid]['lineage_taxid']:
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
                                read.reference_id = refid_dict[read.reference_name]
                                try:
                                        outfile.write(read)
                                except ValueError:
                                    logging.error(f"Error writing read {read.query_name}")
                                    continue
                        except KeyError:
                            # print(read, read.reference_name, read.reference_id)
                            continue
                            
        outfile.close()
    samfile.close()


def write_bam_tags(
    infile,
    outfile,
    total_reads,
    read_taxid_dict,
    taxid_info_dict,
    identity,
    edit_distance,
    minlength,
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
    """
    logging.info(f"* Writing BAM file with taxonomic information tag to {outfile}")
    mode = {"sam": "r", "bam": "rb", "cram": "rc"}
    filetype = infile.split(".")[-1]
    with AlignmentFile(infile, mode[filetype]) as samfile:
        with AlignmentFile(outfile, "wb", template=samfile) as outfile:
            for read in tqdm(samfile, unit="reads", total=total_reads):
                if read.has_tag("NM") and not read.is_unmapped:
                    mismatch = read.get_tag("NM")
                    alnLen = read.query_alignment_length
                    readLen = read.query_length
                    ident = (alnLen - mismatch) / readLen
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
                outfile.write(read)
        outfile.close()
    samfile.close()
