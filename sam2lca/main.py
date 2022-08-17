from sam2lca.parse_alignment import Alignment
from sam2lca.acc2tax import get_mapping
from sam2lca.lca import compute_lca
from sam2lca import utils
from sam2lca.taxonomy import setup_taxopy_db, load_taxonomy_db
from sam2lca.write_alignment import write_bam_tags, write_bam_by_taxid
from pathlib import Path
import logging
from multiprocessing import cpu_count
from tqdm.contrib.concurrent import thread_map
from functools import partial
import sys
import os
from pprint import pprint


def sam2lca(
    sam,
    output=None,
    dbdir=f"{str(Path.home())}/.sam2lca",
    taxonomy="ncbi",
    acc2tax="nucl",
    process=2,
    identity=0.8,
    distance=None,
    length=30,
    conserved=False,
    bam_out=False,
    bam_split_rank=False,
    bam_split_read=50,
):
    """Performs LCA on SAM/BAM/CRAM alignment file

    Args:
        sam (str): Path to SAM/BAM/CRAM alignment file
        output(str): Path to sam2lca output file
        dbdir (str): Path to database storing directory
        taxonomy(str): Type of Taxonomy database
        acc2tax(str): Type of acc2tax database
        process (int): Number of process for parallelization
        identity(float): Minimum alignment identity threshold
        edit_distance(int): Maximum edit distance threshold
        length(int): Minimum alignment length
        bam_out(bool): Write BAM output file with XT tag for TAXID
        bam_split_rank(str): Rank to split BAM output file by TAXID
        bam_split_read(int): Minimum number of reads to split BAM output file by TAXID
    """
    nb_steps = 8 if conserved else 7
    process = cpu_count() if process == 0 else process

    Path(dbdir).mkdir(exist_ok=True)

    avail_taxo, avail_acc2tax = list_available_db(dbdir)
    if taxonomy not in avail_taxo:
        logging.info(
            f"{taxonomy} taxonomy database is not installed. I'm trying to create it for you"
        )
        update_database(taxonomy=taxonomy, dbdir=dbdir, acc2tax=None)
    if acc2tax not in avail_acc2tax:
        logging.info(
            f"{acc2tax} acc2tax database is ls not installed. I'm trying to create it for you"
        )
        update_database(acc2tax=acc2tax, dbdir=dbdir)
    logging.info(f"Step 1/{nb_steps}: Loading taxonomy database")
    TAXDB = load_taxonomy_db(dbdir, taxonomy)

    acc2tax_db = f"{dbdir}/{acc2tax}.db"
    p = Path(acc2tax_db)
    if not p.exists():
        logging.error(
            f"\t'{acc2tax_db}' database seems to be missing, please check your inputs"
        )
        sys.exit(1)
    if not output:
        output = utils.output_file(sam)
    else:
        output = utils.output_file(output)
    utils.check_extension(sam)
    al = Alignment(al_file=sam, nb_steps=nb_steps)
    al.get_refs_taxid(acc2tax_db, TAXDB)
    read_dict = al.get_reads(
        process=process,
        identity=identity,
        edit_distance=distance,
        minlength=length,
        check_conserved=conserved,
    )

    reads_taxid_dict = compute_lca(read_dict, process, nb_steps, taxo_db=TAXDB)
    taxid_counts = utils.count_reads_taxid(reads_taxid_dict)
    taxid_info_dict = utils.taxid_to_lineage(
        taxid_counts,
        output["sam2lca"],
        process=process,
        nb_steps=nb_steps,
        taxo_db=TAXDB,
    )


    if bam_out and bam_split_rank:
        logging.info(f"* BAM files split by TAXID at the {bam_split_rank} level")
        taxids_list = utils.filter_taxid_info_dict(
            taxid_info_dict,
            rank=bam_split_rank, 
            minread=bam_split_read
        ).keys()

        write_bam_by_taxid(
            target_taxids=taxids_list,
            infile=sam,
            outfile_base=output["bam"],
            total_reads = al.total_reads,
            read_taxid_dict=reads_taxid_dict,
            acc2tax_dict=al.acc2tax,
            taxid_info_dict=taxid_info_dict,
            identity=identity,
            edit_distance=distance,
            minlength=length,
            process=process
        )

    if bam_out and not bam_split_rank:
        write_bam_tags(
            infile=sam,
            outfile=output["bam"],
            total_reads=al.total_reads,
            read_taxid_dict=reads_taxid_dict,
            taxid_info_dict=taxid_info_dict,
            identity=identity,
            edit_distance=distance,
            minlength=length,
            process=process
        )
    return taxid_info_dict


def update_database(
    dbdir=f"{str(Path.home())}/.sam2lca",
    taxonomy=None,
    taxo_names=None,
    taxo_nodes=None,
    taxo_merged=None,
    acc2tax="nucl",
    acc2tax_json="https://raw.githubusercontent.com/maxibor/sam2lca/master/data/acc2tax.json",
):
    """Performs LCA on SAM/BAM/CRAM alignment file

    Args:

        dbdir (str): Path to database storing directory
        taxonomy(str): Name of Taxonomy database
        names(str): names.dmp file for taxonomy database. None loads the NCBI taxonomy database
        nodes(str): nodes.dmp file for taxonomy database. None loads the NCBI taxonomy database
        merged(str): merged.dmp file for taxonomy database. None loads the NCBI taxonomy database
        acc2tax(str): Type of acc2tax database
        acc2tax_json(str): Path to acc2tax json file
    """

    
    map_config = utils.get_json(json_path=acc2tax_json)

    if acc2tax is not None:
        logging.info(f"* Downloading/updating acc2tax {acc2tax} database")

        get_mapping(map_config=map_config, maptype=acc2tax, dbdir=dbdir)

    if taxonomy:
        logging.info(f"* Setting up {taxonomy} taxonomy database")
        setup_taxopy_db(
            dbdir=dbdir,
            db_type=taxonomy,
            nodes=taxo_nodes,
            names=taxo_names,
            merged=taxo_merged,
        )


def list_available_db(dbdir, verbose=False):
    """List available taxonomy databases

    Args:
        db_dir(str): Path to sam2lca database directory

    Returns:
        list: List of available taxonomy databases
        list: List of available acc2tax databases
    """
    taxo_db = []
    acc2tax_db = []
    if os.path.exists(dbdir):
        for file in os.listdir(dbdir):
            if file.endswith(".pkl"):
                taxo_db.append(file.split(".")[0])
                continue
            if file.endswith(".db"):
                acc2tax_db.append(file.split(".")[0])
                continue
        if verbose:
            logging.info(f"* Available taxonomy databases: {', '.join(taxo_db)}")
            logging.info(f"* Available acc2tax databases: {', '.join(acc2tax_db)}")
        return taxo_db, acc2tax_db
    else:
        logging.error(f"* {dbdir} does not exist, please create sam2lca database first")
        sys.exit(1)


if __name__ == "__main__":
    sam2lca(
        sam=f"{Path(__file__).parent.resolve()}/../tests/data/aligned.sorted.bam",
        conserved=False,
        tree=f"{Path(__file__).parent.resolve()}/../tests/data/taxonomy/test.tree",
        mappings="test",
        process=0,
        bam_out=True,
    )
