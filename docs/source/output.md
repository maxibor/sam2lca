# Output

sam2lca generates a [JSON](https://www.w3schools.com/python/python_json.asp) file as output, with NCBI Taxonomy IDs as keys.

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