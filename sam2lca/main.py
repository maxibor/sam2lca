from sam2lca.sam_pysam import Alignment
from sam2lca.ncbi_to_taxid import get_mapping
from sam2lca.mapfiles import map_db
from sam2lca.sam_ete import compute_lca_multi
from sam2lca.utils import get_script_dir, count_reads_taxid


def sam2lca(sam, mappings, tree, update, process):
    al = Alignment(al_file=sam, filetype='bam')
    read_dict = al.get_reads()
    get_mapping(mappings, update)
    acc2tax_db = f"{get_script_dir()}/mappings/{map_db[mappings]}"
    reads_taxid_dict = compute_lca_multi(read_dict,
                                         acc2tax_db,
                                         tree,
                                         update,
                                         process)
    print(count_reads_taxid(reads_taxid_dict))
