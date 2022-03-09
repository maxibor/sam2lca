# Databases

The databases available to use with `sam2lca` are:

## Nucleotide databases

- `nucl` for nucleotide/DNA sequences, made of:
    - `nucl_wgs` : nucleotide sequence records of type WGS or TSA
    - `nucl_gb` : nucleotide sequence records that are not WGS or TSA
- `plant_markers` for plant identication based on plant specific markers, made of:
    - `angiosperms353` : Angiosperms353 marker data extracted from [treeoflife.kew.org](https://treeoflife.kew.org/) with sequence headers reformatted as following:

    _Original fasta header_
    ```
    >5821 Gene_Name:dph5 Species:Cyperus_laevigatus Repository:INSDC Sequence_ID:ERR3650073
    ```

    _Reformatted fasta header_
    ```
    >5821_Cyperus_laevigatus Gene_Name:REV7  Repository:INSDC Sequence_ID:ERR3650073
    ```

    > This reformating is necessary to ensure the uniqueness of sequence identifiers. The `fasta` file with reformatted headers (dumped from [treeoflife.kew.org](https://treeoflife.kew.org/) on October 21st, 2021) is available for download here: [angiosperms353_markers.fa.gz](https://edmond.mpdl.mpg.de/api/access/datafile/101862)

    - `ITS` : ITS plant markers data extracted from the [planITS project](https://github.com/apallavicini/PLANiTS). The ITS database is available [ITS.fa.gz](https://edmond.mpdl.mpg.de/api/access/datafile/101863)

    - `rbcL`: rbcL plant marker extraced from [10.3732/apps.1600110](https://doi.org/10.3732/apps.1600110), using the version updated on 09.07.2021, shared by the authors [here](https://figshare.com/collections/rbcL_reference_library_July_2021/5504193). 
    Fasta headers were rewritten to ensure the uniqueness of sequence identifiers and the dabase is available [rbcl.fa.gz](https://edmond.mpdl.mpg.de/api/access/datafile/101864).

    _Original fasta header_

    ```
    >123456 Grabowskia glauca
    ```
    
    __Reformatted fasta header_
    ```
    >rbcL_0_Grabowskia_glauca
    ```

    


## Protein databases

- `prot` for protein sequences, made of:
    - `prot` : protein sequence records which have GI identifiers
    - `pdb` : protein sequence records from the Protein Data Bank

## Test database

- `test` : local database to test sam2lca