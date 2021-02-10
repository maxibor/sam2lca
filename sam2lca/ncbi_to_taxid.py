import pickle
import sys
import hashlib
import os
import urllib.request as urllib
from sam2lca.mapfiles import mapfiles, mapmd5, map_db
from sam2lca.utils import get_script_dir
from sam2lca.config import OPTS_create
from xopen import xopen
from tqdm import tqdm
from subprocess import check_output
import rocksdb
import shutil

class TqdmUpTo(tqdm):
    """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""

    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)  # will also set self.n = b * bsize


def wc(filename, compressed=True):
    if compressed:
        cmd = f"zcat < {filename} | wc -l"
    else:
        cmd = f"wc -l {filename}"
    return int(check_output(cmd, shell=True).split()[0])


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def dl_mappings(mapping_url, md5_url, dbdir):
    """Download mapping file

    Args:
        mapping_url (str): URL of mapping file
        md5_url (str): URL of md5 mapping file
        dbdir (str): path to mapping local directory
    Returns:
        mapping_fname (str): mapping filename 
    """
    testdir = get_script_dir() + '/../tests/data/taxonomy/'

    md5_fname = md5_url.split("/")[-1]
    mapping_fname = mapping_url.split("/")[-1]
    if mapping_fname != 'test.accession2taxid.gz':
        urllib.urlretrieve(md5_url, f"{dbdir}/{md5_fname}")

        with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=mapping_fname) as t:
            urllib.urlretrieve(mapping_url, filename=f"{dbdir}/{mapping_fname}",
                            reporthook=t.update_to, data=None)
    else:
        shutil.copy(testdir+mapping_fname, f"{dbdir}/{mapping_fname}")
        shutil.copy(testdir+md5_fname, f"{dbdir}/{md5_fname}")
    with open(f"{dbdir}/{md5_fname}", 'r') as hf:
        for line in hf:
            md5_hash = line.rstrip().split()[0]
    file_hash = md5(f"{dbdir}/{mapping_fname}")
    try:
        assert file_hash == md5_hash
    except AssertionError:
        print(f"Error while downloading {mapping_url}, please try again")
        sys.exit(1)
    return(mapping_fname)


def mapping_file_to_db(mapdb, mapfile, mapmd5, dbdir):
    """Read mapping to dict and pickle

    Args:
        mapdb (str): Mapping db name
        mapfile(str): Mapping file
        mapmd5 (str): Mapping file md5
        dbdir (str): Mapping file location
    """

    def get_batch_size(nb_entries, nb_batch):
        """Compute of entries to include per batch

        Args:
            nb_entries (int): Number of initial entries
            nb_batch (int): Desired number of batches
        Returns:
            (int): Batch size
        """
        batch_size = min(nb_entries, int(nb_entries/nb_batch))
        if batch_size == 0:
            return(nb_entries)
        else:
            return(batch_size)

    batch = rocksdb.WriteBatch()

    nlines = wc(f"{dbdir}/{mapfile}")
    nb_batch = 100
    
    batch_size = get_batch_size(nlines-1, nb_batch)

    with xopen(f"{dbdir}/{mapfile}") as acc:
        i = 0
        for line in tqdm(acc, total=nlines):
            if not line.startswith("accession"):
                acc_line = line.split()
                acc_num = acc_line[1]
                acc_taxid = acc_line[2]
                batch.put(bytes(acc_num, encoding='utf8'),
                        bytes(acc_taxid, encoding='utf8'))
            if i % batch_size == 0:
                DATABASE.write(batch)
                batch = rocksdb.WriteBatch()
    DATABASE.write(batch)
    os.remove(f"{dbdir}/{mapfile}")
    os.remove(f"{dbdir}/{mapmd5}")


def get_mapping(maptype, update, dbdir):
    """Get genbank to taxid mapping

    Args:
        maptype(str): Mapping type: nucl or prot
    Returns:
        dict: 
    """
    if maptype not in mapfiles.keys():
        print(f"mapping type not in {', '.join(list(mapfiles.keys()))}")
        sys.exit(1)

    os.makedirs(dbdir, exist_ok=True)

    if map_db[maptype] not in os.listdir(dbdir) or update:
        global DATABASE
        db_name = f"{dbdir}/{map_db[maptype]}"
        DATABASE = rocksdb.DB(db_name, OPTS_create)

        nb_files = len(mapfiles[maptype])
        for i in range(nb_files):
            mapfile = mapfiles[maptype][i].split("/")[-1]
            if mapfile not in os.listdir(dbdir) or update:
                try:
                    os.remove(f"{dbdir}/{mapfile}")
                except FileNotFoundError as e:
                    print(f"No local copy of {mapfile}")

                print(f"Downloading '{maptype}' acc2tax mapping file [{i+1}/{nb_files}]")
                mapname = dl_mappings(mapping_url=mapfiles[maptype][i],
                                    md5_url=mapmd5[maptype][i], dbdir=dbdir)
            else:
                mapname = mapfiles[maptype][i].split("/")[-1]
            print("Inserting mappings into database")
            mapping_file_to_db(mapfile=mapname,
                            mapdb=map_db[maptype],
                            dbdir=dbdir,
                            mapmd5=mapmd5[maptype][i].split("/")[-1])

    # DB = rocksdb.DB(f"{dbdir}/{map_db[maptype]}", opts=rocksdb.Options(), read_only=True)
    # return(f"{dbdir}/{map_db[maptype]}")


if __name__ == "__main__":
    allmaps = get_mapping('pdb')
