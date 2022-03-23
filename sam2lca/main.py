from sam2lca.sam_pysam import Alignment
from sam2lca.ncbi_to_taxid import get_mapping
from sam2lca.mapfiles import map_db
from sam2lca.sam_ete import compute_lca_multi
from sam2lca import utils
from sam2lca.config import NCBI
from sam2lca.write_alignment import write_bam_tags
from pathlib import Path
import logging
from multiprocessing import cpu_count


def sam2lca(
    sam,
    mappings,
    output=None,
    tree=None,
    dbdir=f"{str(Path.home())}/.sam2lca",
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
        output(str): Path to sam2lca output file
        tree (str): Path to optional taxonomic tree
        dbdir (str): Path to database storing directory
        process (int): Number of process for parallelization
        identity(float): Minimum identity
        length(int): Minimum alignment length
        bam_out(bool): Write BAM output file with XT tag for TAXID
    """
    nb_steps = 8 if conserved else 7
    process = cpu_count() if process == 0 else process

    acc2tax_db = f"{dbdir}/{map_db[mappings]}"
    p = Path(acc2tax_db)
    if not p.exists():
        logging.info(
            f"\t'{acc2tax_db}' database seems to be missing, I'm creating it for you."
        )
        get_mapping(mappings, dbdir=dbdir, update=False)

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

    reads_taxid_dict = compute_lca_multi(read_dict, tree, process, nb_steps)
    taxid_counts = utils.count_reads_taxid(reads_taxid_dict)
    utils.taxid_to_lineage(
        taxid_counts, output["sam2lca"], process=process, nb_steps=nb_steps
    )
    if bam_out:
        write_bam_tags(sam, output["bam"], reads_taxid_dict)


def update_database(mappings, dbdir, ncbi):
    """Performs LCA on SAM/BAM/CRAM alignment file

    Args:
        mappings (str): Type of Acc2Tax mapping
        dbdir (str): Path to database stroring directory
        ncbi (bool): Updates NCBI taxonomic tree
    """
    logging.info("Downloading/updating database")
    NCBI
    if ncbi:
        NCBI.update_taxonomy_database()
    get_mapping(mappings, dbdir=dbdir, update=True)


if __name__ == "__main__":
    sam2lca(
        sam=f"{Path(__file__).parent.resolve()}/../tests/data/aligned.sorted.bam",
        conserved=False,
        tree=f"{Path(__file__).parent.resolve()}/../tests/data/taxonomy/test.tree",
        mappings="test",
        process=0,
        bam_out=True,
    )
