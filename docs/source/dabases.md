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

    > This reformating is necessary to ensure the uniqueness of sequence identifiers. The `fasta` file with reformatted headers (dumped from [treeoflife.kew.org](https://treeoflife.kew.org/) on October 21st, 2021) is available [for download here](https://edmond.mpdl.mpg.de/imeji/file/ee/02/37/69-b50d-4f5b-8729-efacb3ac7c13/0/original/b56871efbefaf0ddad48ce12a59a4c54.gz?filename=angiosperms353_markers.fa.gz).

    - `ITS` : ITS plant markers data extracted from the [planITS project](https://github.com/apallavicini/PLANiTS). The ITS database is available [for download here](https://edmond.mpdl.mpg.de/imeji/file/75/07/fa/81-1055-4439-977a-647584734e91/0/original/03d581a74d49d3ebe044467a1350220c.gz?filename=ITS.fa.gz)

    - `rbcL`: rbcL plant marker extraced from [10.3732/apps.1600110](https://doi.org/10.3732/apps.1600110), using the version updated on 09.07.2021, shared by the authors [here](https://figshare.com/collections/rbcL_reference_library_July_2021/5504193). 
    Fasta headers were rewritten to ensure the uniqueness of sequence identifiers and the dabase is available [here](https://edmond.mpdl.mpg.de/imeji/file/5c/4d/91/3a-c28e-4cbc-82b6-b0bf80f1cda9/0/original/3a8d6c5884c45c30f2435e055861929e.gz?filename=rbcl.fa.gz)

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