import pkg_resources

mapfiles = {
    'nucl': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz'],
    'prot': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz'],
    'test' : [pkg_resources.resource_filename('sam2lca.data', 'test.accession2taxid.gz')]
}

mapmd5 = {
    'nucl': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5'],
    'prot': ['https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz.md5', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5'],
    'test' : [pkg_resources.resource_filename('sam2lca.data', 'test.accession2taxid.gz.md5')]
}

map_db = {
    'nucl': 'nucl.db',
    'prot': 'prot.db',
    'test': 'test.db'
}
