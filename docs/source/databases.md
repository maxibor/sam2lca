# Databases

sam2lca uses two different type of databases:

- a **taxonomy** database to infer the Lowest Common Ancestor (LCA) and retrieve the names and lineage associated to a *taxonomic identifier* (TAXID)
- an **acc2tax** or *accession to TAXID* database, to match sequence accession to a taxonomic identifier

For each of these databases, sam2lca offers different possibilities.

## Taxonomy databases

## `ncbi`

If you're not sure what to use, stick with `ncbi`

## `gtdb`

Not yet natively implemented, but see the `custom` taxonomy database below and have a look at [gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump)

## `custom`

You can provide your own taxonomy database by providing the following files

- `names.dmp`
- `nodes.dmp`
- `merged.dmp`

## acc2tax - *accession to TAXID* databases

### Nucleotide databases

- `nucl` for nucleotide/DNA sequences, made of:
  - *nucl_wgs* : nucleotide sequence records of type WGS or TSA
  - *nucl_gb* : nucleotide sequence records that are not WGS or TSA
- `plant_markers` for plant identication based on plant specific markers, made of:
  - *angiosperms353* : Angiosperms353 marker data extracted from [treeoflife.kew.org](https://treeoflife.kew.org/) with sequence headers reformatted as following:

    *Original fasta header*

    ```
    >5821 Gene_Name:dph5 Species:Cyperus_laevigatus Repository:INSDC Sequence_ID:ERR3650073
    ```

    *Reformatted fasta header*

    ```
    >5821_Cyperus_laevigatus Gene_Name:REV7  Repository:INSDC Sequence_ID:ERR3650073
    ```

    > This reformating is necessary to ensure the uniqueness of sequence identifiers. The `fasta` file with reformatted headers (dumped from [treeoflife.kew.org](https://treeoflife.kew.org/) on October 21st, 2021) is available for download here: [angiosperms353_markers.fa.gz](https://edmond.mpdl.mpg.de/api/access/datafile/101862)

  - *ITS* : ITS plant markers data extracted from the [planITS project](https://github.com/apallavicini/PLANiTS). The ITS database is available [ITS.fa.gz](https://edmond.mpdl.mpg.de/api/access/datafile/101863)

  - *rbcL*: rbcL plant marker extraced from [10.3732/apps.1600110](https://doi.org/10.3732/apps.1600110), using the version updated on 09.07.2021, shared by the authors [here](https://figshare.com/collections/rbcL_reference_library_July_2021/5504193).
    Fasta headers were rewritten to ensure the uniqueness of sequence identifiers and the dabase is available [rbcl.fa.gz](https://edmond.mpdl.mpg.de/api/access/datafile/101864).

    *Original fasta header*

    ```
    >123456 Grabowskia glauca
    ```

    _*Reformatted fasta header*

    ```
    >rbcL_0_Grabowskia_glauca
    ```

### Protein databases

- `prot` for protein sequences, made of:
  - *prot* : protein sequence records which have GI identifiers
  - *pdb* : protein sequence records from the Protein Data Bank

### Test database

- `test` : local database to test sam2lca

### Custom database

With sam2lca, you can provide a custom database to map accession numbers to TAXIDs.

To do so, sam2lca can accept a [`JSON`](https://www.json.org/json-en.html) file, with the `--acc2tax_json` flag in the `sam2lca update-db` subcommand in combination with `--acc2tax custom`.

For example:

```bash
sam2lca update-db --acc2tax custom --acc2tax_json acc2tax.json
```

This JSON file should be formatted as below:

```json
{
    "mapfiles": {
        "[name_of_mapping]": [
            "path/url_to_compressed_accession2taxid.gz file"
        ]
    },
    "mapmd5": {
        "[name_of_mapping]": [
            "path/url_to_compressed_accession2taxid.gz md5sumfile"
        ]
    },
    "map_db": {
        "[name_of_mapping]": "Name of custom.db"
    }
}
```

An example `json` file  can be found here: [map_config.json](https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/acc2tax.json)
