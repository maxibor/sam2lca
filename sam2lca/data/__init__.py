import pkg_resources
import json

acc2tax_default = {
    "mapfiles": {
        "nucl": [
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
        ],
        "prot": [
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz",
        ],
        "plant_markers": [
            "https://raw.githubusercontent.com/maxibor/sam2lca/master/sam2lca/data/angiosperms353.accession2taxid.gz",
            "https://raw.githubusercontent.com/maxibor/sam2lca/master/sam2lca/data/ITS.accession2taxid.gz",
            "https://raw.githubusercontent.com/maxibor/sam2lca/master/sam2lca/data/rbcl.accession2taxid.gz",
            "https://raw.githubusercontent.com/maxibor/sam2lca/taxopy/sam2lca/data/SILVA_138.1_SSURef_NR99.accession2taxid.gz",
        ],
        "gtdb_r207": [
            "https://raw.githubusercontent.com/maxibor/sam2lca/taxopy/sam2lca/data/gtdb_r207.accession2taxid.gz"
        ],
        "test": [pkg_resources.resource_filename(__name__, "test.accession2taxid.gz")],
    },
    "mapmd5": {
        "nucl": [
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5",
        ],
        "prot": [
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz.md5",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5",
        ],
        "plant_markers": [
            "https://raw.githubusercontent.com/maxibor/sam2lca/master/sam2lca/data/angiosperms353.accession2taxid.gz.md5",
            "https://raw.githubusercontent.com/maxibor/sam2lca/master/sam2lca/data/ITS.accession2taxid.gz.md5",
            "https://raw.githubusercontent.com/maxibor/sam2lca/master/sam2lca/data/rbcl.accession2taxid.gz.md5",
            "https://raw.githubusercontent.com/maxibor/sam2lca/taxopy/sam2lca/data/SILVA_138.1_SSURef_NR99.accession2taxid.gz.md5",
        ],
        "gtdb_r207": [
            "https://raw.githubusercontent.com/maxibor/sam2lca/taxopy/sam2lca/data/gtdb_r207.accession2taxid.gz.md5"
        ],
        "test": [
            pkg_resources.resource_filename(__name__, "test.accession2taxid.gz.md5")
        ],
    },
    "map_db": {
        "nucl": "nucl.db",
        "prot": "prot.db",
        "plant_markers": "plant_markers.db",
        "gtdb_r207": "gtdb_r207.db",
        "test": "test.db",
    },
}


def get_map_config(
    map_config_file,
    base_config=acc2tax_default,
):
    """Reads map_config.json file
    Args:
        map_config_file (str): Path to acc2tax.json file
        base_config (dict): Default map_config
    Returns:
        dict: map_config
        acc2tax_name: name of new acc2tax db
    """

    with open(map_config_file, "r") as f:
        config = json.load(f)
    acc2tax_name = list(config["mapfiles"].keys())[0]
    for k in base_config:
        base_config[k].update(config[k])

    return base_config, acc2tax_name
