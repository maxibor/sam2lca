from sam2lca.utils import get_script_dir

mapfiles = {
    'nucl': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz'],
    'prot': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz'],
    'test': [get_script_dir()+'/../tests/data/taxonomy/test.accession2taxid.gz']
}

mapmd5 = {
    'nucl': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5'],
    'prot': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz.md5', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5'],
    'test': [get_script_dir()+'/../tests/data/taxonomy/test.accession2taxid.gz.md5']
}

map_db = {
    'nucl': 'nucl.db',
    'prot': 'prot.db',
    'test': 'test.db'
}
