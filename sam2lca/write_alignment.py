import logging
from pysam import AlignmentFile


def write_bam_tags(infile, outfile, read_taxid_dict):
    """Write alignment file with taxonomic information tag

    Args:
        infile (str): Path to alignment file
        outfile (str): Path to output file
        read_taxid_dict (dict): Dictionary of read with their corresponding LCA taxid
    """
    logging.info(f"Step 7/6: Writing BAM file with XT TAXID tags to {outfile}")
    mode = {"sam": "r", "bam": "rb", "cram": "rc"}
    filetype = infile.split(".")[-1]
    with AlignmentFile(infile, mode[filetype]) as samfile:
        with AlignmentFile(outfile, "wb", template=samfile) as outfile:
            for read in samfile:
                read.set_tag("XT", read_taxid_dict[read.query_name], "i")
                outfile.write(read)
        outfile.close()
    samfile.close()
