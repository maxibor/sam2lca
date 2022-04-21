import sys
import os
from sam2lca.rocksdb import OPTS_create
from sam2lca.utils import download_and_checksum, wc
from xopen import xopen
from tqdm import tqdm
import rocksdb
import logging


def mapping_file_to_db(mapfile, mapmd5, dbdir):
    """Read mapping to dict and pickle

    Args:
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
        batch_size = min(nb_entries, int(nb_entries / nb_batch))
        if batch_size == 0:
            return nb_entries
        else:
            return batch_size

    batch = rocksdb.WriteBatch()

    nlines = wc(f"{dbdir}/{mapfile}")
    nb_batch = 100

    batch_size = get_batch_size(nlines - 1, nb_batch)

    with xopen(f"{dbdir}/{mapfile}") as acc:
        i = 0
        for line in tqdm(acc, total=nlines):
            if not line.startswith("accession"):
                acc_line = line.split()
                acc_num = acc_line[1]
                acc_taxid = acc_line[2]
                batch.put(
                    bytes(acc_num, encoding="utf8"), bytes(acc_taxid, encoding="utf8")
                )
            if i % batch_size == 0:
                DATABASE.write(batch)
                batch = rocksdb.WriteBatch()
    DATABASE.write(batch)
    os.remove(f"{dbdir}/{mapfile}")
    os.remove(f"{dbdir}/{mapmd5}")


def get_mapping(map_config, maptype, dbdir):
    """Get genbank to taxid mapping

    Args:
        map_config (dict): Mapping config
        maptype(str): Mapping type
        dbdir(str): Directory to store taxonomy databases
    """
    mapfiles = map_config["mapfiles"]
    mapmd5 = map_config["mapmd5"]
    map_db = map_config["map_db"]

    if maptype not in mapfiles.keys():
        logging.error(f"mapping type not in {', '.join(list(mapfiles.keys()))}")
        sys.exit(1)

    os.makedirs(dbdir, exist_ok=True)

    global DATABASE
    db_name = f"{dbdir}/{map_db[maptype]}"
    DATABASE = rocksdb.DB(db_name, OPTS_create)

    nb_files = len(mapfiles[maptype])
    for i in range(nb_files):
        mapfile = mapfiles[maptype][i].split("/")[-1]
        try:
            os.remove(f"{dbdir}/{mapfile}")
        except FileNotFoundError as e:
            pass
        logging.info(f"Downloading '{maptype}' acc2tax mapping file [{i+1}/{nb_files}]")
        mapname = download_and_checksum(
            file_url=mapfiles[maptype][i],
            md5_url=mapmd5[maptype][i],
            dbdir=dbdir,
        )
        logging.info("Inserting mappings into database")
        mapping_file_to_db(
            mapfile=mapname,
            dbdir=dbdir,
            mapmd5=mapmd5[maptype][i].split("/")[-1],
        )

    # DB = rocksdb.DB(f"{dbdir}/{map_db[maptype]}", opts=rocksdb.Options(), read_only=True)
    # return(f"{dbdir}/{map_db[maptype]}")


if __name__ == "__main__":
    allmaps = get_mapping("pdb")
