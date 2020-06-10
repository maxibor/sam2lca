import pickle
import sys
import hashlib
import os
import urllib.request as urllib
from sam2lca.mapfiles import mapfiles, mapmd5, map_db
from sam2lca.utils import get_script_dir
from xopen import xopen
from tqdm import tqdm
from subprocess import check_output
import rocksdb


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


def dl_mappings(mapping_url, md5_url, mapdir):
    """Download mapping file

    Args:
        mapping_url (str): URL of mapping file
        md5_url (str): URL of md5 mapping file
        mapdir (str): path to mapping local directory
    Returns:
        mapping_fname (str): mapping filename 
    """
    md5_fname = md5_url.split("/")[-1]
    mapping_fname = mapping_url.split("/")[-1]
    urllib.urlretrieve(md5_url, f"{mapdir}/{md5_fname}")

    with TqdmUpTo(unit='B', unit_scale=True, miniters=1, desc=mapping_fname) as t:
        urllib.urlretrieve(mapping_url, filename=f"{mapdir}/{mapping_fname}",
                           reporthook=t.update_to, data=None)

    with open(f"{mapdir}/{md5_fname}", 'r') as hf:
        for line in hf:
            md5_hash = line.rstrip().split()[0]
    file_hash = md5(f"{mapdir}/{mapping_fname}")
    try:
        assert file_hash == md5_hash
    except AssertionError:
        print(f"Error while downloading {mapping_url}, please try again")
        sys.exit(1)
    return(mapping_fname)


def mapping_file_to_db(mapdb, mapfile, mapdir):
    """Read mapping to dict and pickle

    Args:
        mapdb (str): Mapping db name
        mapfile(str): Mapping file
        mapdir (str): Mapping file location
    """
    db_name = f"{mapdir}/{mapdb}"
    db = rocksdb.DB(db_name, rocksdb.Options(create_if_missing=True))
    batch = rocksdb.WriteBatch()

    nlines = wc(f"{mapdir}/{mapfile}")
    with xopen(f"{mapdir}/{mapfile}") as acc:
        for line in tqdm(acc, total=nlines):
            if not line.startswith("accession"):
                acc_line = line.split()
                acc_num = acc_line[1]
                acc_taxid = acc_line[2]
                batch.put(bytes(acc_num, encoding='utf8'),
                          bytes(acc_taxid, encoding='utf8'))

    db.write(batch)


def get_mapping(maptype, update):
    """Get genbank to taxid mapping

    Args:
        maptype(str): Mapping type: nucl or prot
    Returns:
        dict: 
    """
    if maptype not in mapfiles.keys():
        print(f"mapping type not in {', '.join(list(mapfiles.keys()))}")
        sys.exit(1)

    mapdir = get_script_dir()+"/mappings"

    if map_db[maptype] not in os.listdir(mapdir) or update:
        mapfile = mapfiles[maptype].split("/")[-1]
        if mapfile not in os.listdir(mapdir) or update:
            try:
                os.remove(f"{mapdir}/{mapfile}")
            except FileNotFoundError as e:
                print(f"No local copy of {mapfile}")

            print(f"Downloading {maptype}  acc2tax mapping file")
            mapname = dl_mappings(mapping_url=mapfiles[maptype],
                                  md5_url=mapmd5[maptype], mapdir=mapdir)
        else:
            mapname = mapfiles[maptype].split("/")[-1]
        print("Inserting mappings into database")
        mapping_file_to_db(mapfile=mapname,
                           mapdb=map_db[maptype],
                           mapdir=mapdir)
    # else:
    #     print(f"{maptype} mapping file present, reading from pickle file")
    #     with open(f"{mapdir}/{map_pickle[maptype]}", 'rb') as p:
    #         mapdict = pickle.load(p)
    # return(mapdict)


if __name__ == "__main__":
    allmaps = get_mapping('pdb')
