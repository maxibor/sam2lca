import logging
from pysam import AlignmentFile


def write_bam_tags(infile, outfile, read_taxid_dict, identity, minlength, unclassified_taxid=12908):
    """Write alignment file with taxonomic information tag

    Args:
        infile (str): Path to alignment file
        outfile (str): Path to output file
        read_taxid_dict (dict): Dictionary of read with their corresponding LCA taxid
        identity (float): Minimum identity for alignment
        minlength (int): Minimum length for alignment
    """
    logging.info(f"* BAM file, with XT tags set to LCA TAXID, to {outfile}")
    mode = {"sam": "r", "bam": "rb", "cram": "rc"}
    filetype = infile.split(".")[-1]
    with AlignmentFile(infile, mode[filetype]) as samfile:
        with AlignmentFile(outfile, "wb", template=samfile) as outfile:
            for read in samfile:
                if read.has_tag("NM") and not read.is_unmapped:
                    mismatch = read.get_tag("NM")
                    alnLen = read.query_alignment_length
                    readLen = read.query_length
                    ident = (alnLen - mismatch) / readLen
                    if ident >= identity and alnLen >= minlength:
                        try:
                            read.set_tag("XT", read_taxid_dict[read.query_name], "i")
                        except KeyError:
                            read.set_tag("XT", unclassified_taxid, "i")
                    else:
                        read.set_tag("XT", unclassified_taxid, "i")
                else:
                    read.set_tag("XT", unclassified_taxid, "i")
                outfile.write(read)
        outfile.close()
    samfile.close()
