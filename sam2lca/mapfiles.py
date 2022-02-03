import pkg_resources

mapfiles = {
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
    ],
    "test": [pkg_resources.resource_filename(__name__, "data/test.accession2taxid.gz")],
}

mapmd5 = {
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
    ],
    "test": [
        pkg_resources.resource_filename(__name__, "data/test.accession2taxid.gz.md5")
    ],
}

map_db = {
    "nucl": "nucl.db",
    "prot": "prot.db",
    "plant_markers": "plant_markers.db",
    "test": "test.db",
}
