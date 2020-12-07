from sam2lca.sam_pysam import Alignment
from sam2lca.ncbi_to_taxid import get_mapping
from sam2lca.mapfiles import map_db
from sam2lca.sam_ete import compute_lca_multi
from sam2lca import utils
from sam2lca.config import NCBI

def sam2lca(sam, mappings, tree, update, process, identity, length, dbdir, output):
    """Performs LCA on SAM/BAM/CRAM alignment file

    Args:
        sam (str): Path to SAM/BAM/CRAM alignment file
        mappings (str): Type of Acc2Tax mapping
        tree (str): Optional taxonomic tree
        update (bool): Update NCBI taxonomy and Acc2Tax
        process (int): Number of process for parallelization
        identity(float): Minimum identity
        length(int): Minimum alignment length
        output(str): Path to sam2lca output file
    """
    if update:
        NCBI.update_taxonomy_database()
    if not output:
        output = utils.output_file(sam)
    utils.check_extension(sam)
    al = Alignment(al_file=sam)
    read_dict = al.get_reads(
        process=process, identity=identity, minlength=length)
    get_mapping(mappings, update, dbdir=dbdir)
    acc2tax_db = f"{dbdir}/{map_db[mappings]}"
    reads_taxid_dict = compute_lca_multi(read_dict,
                                         acc2tax_db,
                                         tree,
                                         update,
                                         process)
    taxid_counts = utils.count_reads_taxid(reads_taxid_dict)
    utils.taxid_to_lineage(taxid_counts, output)
