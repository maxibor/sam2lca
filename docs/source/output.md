# Output

sam2lca generates:

- a [`JSON`](https://www.w3schools.com/python/python_json.asp) outfile
- a `CSV` output file
- (optionally), a `BAM` alignment file with the `XT` tag set to the NCBI Taxonomy IDs computed by the LCA.

## JSON

A JSON file with NCBI Taxonomy IDs as keys.

- **name**: scientific name of the taxon
- **rank**: taxonomic rank of the taxon
- **count**: number of reads mapping to the taxon
- **lineage**: taxonomic lineage of the taxon

Example:

```bash
{
    "543": {
        "name": "Enterobacteriaceae", 
        "rank": "family", 
        "count": 2152, 
        "lineage": [
            {
                "no rank": "root"
            }, 
            {
                "no rank": "cellular organisms"
            }, 
            {
                "superkingdom": "Bacteria"
            }, 
            {
                "phylum": "Proteobacteria"
            }, 
            {
                "class": "Gammaproteobacteria"
            }, 
            {
                "order": "Enterobacterales"
            }, 
            {
                "family": "Enterobacteriaceae"
            }
        ]
    }, 
    "300267": {
        "name": "Shigella dysenteriae Sd197", 
        "rank": "no rank", 
        "count": 338, 
        "lineage": [
            {
                "no rank": "root"
            }, 
            {
                "no rank": "cellular organisms"
            }, 
            {
                "superkingdom": "Bacteria"
            }, 
            {
                "phylum": "Proteobacteria"
            }, 
            {
                "class": "Gammaproteobacteria"
            }, 
            {
                "order": "Enterobacterales"
            }, 
            {
                "family": "Enterobacteriaceae"
            }, 
            {
                "genus": "Shigella"
            }, 
            {
                "species": "Shigella dysenteriae"
            }, 
            {
                "no rank": "Shigella dysenteriae Sd197"
            }
        ]
    }, 
    "511145": {
        "name": "Escherichia coli str. K-12 substr. MG1655", 
        "rank": "no rank", 
        "count": 385, 
        "lineage": [
            {
                "no rank": "root"
            }, 
            {
                "no rank": "cellular organisms"
            }, 
            {
                "superkingdom": "Bacteria"
            }, 
            {
                "phylum": "Proteobacteria"
            }, 
            {
                "class": "Gammaproteobacteria"
            }, 
            {
                "order": "Enterobacterales"
            }, 
            {
                "family": "Enterobacteriaceae"
            }, 
            {
                "genus": "Escherichia"
            }, 
            {
                "species": "Escherichia coli"
            }, 
            {
                "no rank": "Escherichia coli K-12"
            }, 
            {
                "no rank": "Escherichia coli str. K-12 substr. MG1655"
            }
        ]
    }
}
```

## CSV

**Rows**: Taxons

**Columns**:

- `TAXID`: NCBI taxonomy ID
- `name`: Name of the taxon
- `rank`: Taxonomic rank
- `count`: Number of reads assigned to this taxon
- `lineage`: Taxonomic lineage of this taxon, each taxonomic level being separated by a `-` sign.

```python
TAXID,name,rank,count,lineage
543,Enterobacteriaceae,family,2152,'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Bacteria' - 'phylum': 'Proteobacteria' - 'class': 'Gammaproteobacteria' - 'order': 'Enterobacterales' - 'family': 'Enterobacteriaceae'
511145,Escherichia coli str. K-12 substr. MG1655,no rank,385,'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Bacteria' - 'phylum': 'Proteobacteria' - 'class': 'Gammaproteobacteria' - 'order': 'Enterobacterales' - 'family': 'Enterobacteriaceae' - 'genus': 'Escherichia' - 'species': 'Escherichia coli' - 'strain': 'Escherichia coli K-12' - 'no rank': 'Escherichia coli str. K-12 substr. MG1655'
300267,Shigella dysenteriae Sd197,strain,338,'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Bacteria' - 'phylum': 'Proteobacteria' - 'class': 'Gammaproteobacteria' - 'order': 'Enterobacterales' - 'family': 'Enterobacteriaceae' - 'genus': 'Shigella' - 'species': 'Shigella dysenteriae' - 'strain': 'Shigella dysenteriae Sd197'
```

## BAM

> Only  generated when running `sam2lca analyze` with the `-b`/`--bam_out` flag

The input alignment file is written as a `bam` file, with an extra tag `XT` (of type `int`/`i`) set to the NCBI Taxonomy IDs attributed to each sequencing read by the LCA algorithm.

```bam
shigella_dysenteriae_1462 355 NC_007606.1 276532 1 69M = 276570 113 AGTGCTTTTGCCGTTACGCACCACCCCGTCAGTAGCCGAACAGGAGGGACAGCTGATAGAAACAGAAGC @+:C22<1:??D:FFC;.6@@A;C;EE)=CEAH?@B(6=C>6;5;;;;><CC:;@BBCCACC>38@C@B AS:i:-2 XS:i:0 XN:i:0 XM:i:1 XO:i:0 XG:i:0 NM:i:1 MD:Z:36T32 YS:i:0 YT:Z:CP XT:i:543
shigella_dysenteriae_61 339 NC_007606.1 15921 6 76M = 15869 -128 CTGAACGCGAGAAGCAGATTGTTTTACGAGCCAAGCGCTTAATGCGGGTGCGCAGCGTCAGGTTATTGCGTTCAAT @>CBA?<CCCACCCCCCCCCAA?C?BBB?CCCCABD?AEEFIIIIIIIGBBHD<IIGGGIIIIIIIIIHGHHHFFD AS:i:-5 XS:i:-5 XN:i:0 XM:i:1 XO:i:0 XG:i:0 NM:i:1 MD:Z:70C5YS:i:0 YT:Z:CP XT:i:300267
```

Reads belonging to a specific TAXID can then be filtered with [`samtools view`](http://www.htslib.org/doc/samtools-view.html) like this: `samtools view --tag XT:[YOUR_TAXID_OF_CHOICE] [YOURFILE.bam]`

Example: `samtools view --tag XT:300267 aligned.sorted.bam`
