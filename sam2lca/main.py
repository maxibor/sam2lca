from glob import glob
from sam2lca.sam_pysam import Alignment
from sam2lca.ncbi_to_taxid import get_mapping
from sam2lca.mapfiles import get_map_config, base_map_config
from sam2lca.lca import compute_lca
from sam2lca import utils
from sam2lca.config import setup_taxopy_db
from sam2lca.write_alignment import write_bam_tags
from pathlib import Path
import logging
from multiprocessing import cpu_count


def sam2lca(
    sam,
    mappings,
    map_config=None,
    output=None,
    dbdir=f"{str(Path.home())}/.sam2lca",
    names=None,
    nodes=None,
    merged=None,
    process=2,
    identity=0.8,
    length=30,
    conserved=False,
    bam_out=False,
):
    """Performs LCA on SAM/BAM/CRAM alignment file

    Args:
        sam (str): Path to SAM/BAM/CRAM alignment file
        mappings (str): Type of Acc2Tax mapping
        map_config(str): Path to map_config json file
        output(str): Path to sam2lca output file
        dbdir (str): Path to database storing directory
        names(str): names.dmp file for taxonomy database. None loads the NCBI taxonomy database
        nodes(str): nodes.dmp file for taxonomy database. None loads the NCBI taxonomy database
        merged(str): merged.dmp file for taxonomy database. None loads the NCBI taxonomy database
        process (int): Number of process for parallelization
        identity(float): Minimum identity
        length(int): Minimum alignment length
        bam_out(bool): Write BAM output file with XT tag for TAXID
    """
    nb_steps = 8 if conserved else 7
    process = cpu_count() if process == 0 else process

    print(names, nodes, merged)

    logging.info(f"Step 1/{nb_steps}: Loading taxonomy database")

    TAXDB = setup_taxopy_db(db_path=dbdir, names=names, nodes=nodes, merged=merged)

    if map_config is None:
        map_config = base_map_config
    else:
        map_config = get_map_config(map_config_file=map_config)

    acc2tax_db = f"{dbdir}/{map_config['map_db'][mappings]}"
    p = Path(acc2tax_db)
    if not p.exists():
        logging.info(
            f"\t'{acc2tax_db}' database seems to be missing, I'm creating it for you."
        )
        get_mapping(map_config=map_config, maptype=mappings, dbdir=dbdir, update=False)

    if not output:
        output = utils.output_file(sam)
    else:
        output = utils.output_file(output)
    utils.check_extension(sam)
    al = Alignment(al_file=sam, nb_steps=nb_steps)
    al.get_refs_taxid(acc2tax_db)
    read_dict = al.get_reads(
        process=process,
        identity=identity,
        minlength=length,
        check_conserved=conserved,
    )

    reads_taxid_dict = compute_lca(read_dict, process, nb_steps, taxo_db=TAXDB)
    taxid_counts = utils.count_reads_taxid(reads_taxid_dict)
    utils.taxid_to_lineage(
        taxid_counts,
        output["sam2lca"],
        process=process,
        nb_steps=nb_steps,
        taxo_db=TAXDB,
    )
    if bam_out:
        write_bam_tags(sam, output["bam"], reads_taxid_dict)


def update_database(mappings, dbdir, update, map_config=None):
    """Performs LCA on SAM/BAM/CRAM alignment file

    Args:
        mappings (str): Type of Acc2Tax mapping
        dbdir (str): Path to database stroring directory
        update (bool): Updates taxonomic database
        map_config(str): Path to map_config json file
    """
    logging.info("* Downloading/updating acc2tax database ")

    if map_config is None:
        map_config = base_map_config
    else:
        map_config = get_map_config(map_config_file=map_config)

    get_mapping(map_config=map_config, maptype=mappings, dbdir=dbdir, update=True)

    logging.info("* Setting up taxonomy database")
    print(dbdir)
    TAXDB = setup_taxopy_db(db_path=dbdir)


if __name__ == "__main__":
    sam2lca(
        sam=f"{Path(__file__).parent.resolve()}/../tests/data/aligned.sorted.bam",
        conserved=False,
        tree=f"{Path(__file__).parent.resolve()}/../tests/data/taxonomy/test.tree",
        mappings="test",
        process=0,
        bam_out=True,
    )
