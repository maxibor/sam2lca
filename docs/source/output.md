# Output

sam2lca generates:

- a [`JSON`](https://www.w3schools.com/python/python_json.asp) file
- a `CSV` file
- (optionally), a `BAM` alignment file with the `XT` tag set to the NCBI Taxonomy IDs computed by the LCA.

## JSON

A JSON file with NCBI Taxonomy IDs as keys.

- `name`: scientific name of the taxon
- `rank`: taxonomic rank of the taxon
- `count_taxon`: number of reads mapping to the taxon
- `count_descendant`: total number of reads belonging to the descendants of the taxon
- `lineage`: taxonomic lineage of the taxon

Example:

```bash
{
{
    "1": {
        "name": "root",
        "rank": "no rank",
        "count_taxon": 0,
        "count_descendant": 2875,
        "lineage": {}
    },
    "2": {
        "name": "Bacteria",
        "rank": "superkingdom",
        "count_taxon": 0,
        "count_descendant": 2875,
        "lineage": {
            "superkingdom": "Bacteria"
        }
    },
    "543": {
        "name": "Enterobacteriaceae",
        "rank": "family",
        "count_taxon": 2152,
        "count_descendant": 2875,
        "lineage": {
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "561": {
        "name": "Escherichia",
        "rank": "genus",
        "count_taxon": 0,
        "count_descendant": 385,
        "lineage": {
            "genus": "Escherichia",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "562": {
        "name": "Escherichia coli",
        "rank": "species",
        "count_taxon": 0,
        "count_descendant": 385,
        "lineage": {
            "species": "Escherichia coli",
            "genus": "Escherichia",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "620": {
        "name": "Shigella",
        "rank": "genus",
        "count_taxon": 0,
        "count_descendant": 338,
        "lineage": {
            "genus": "Shigella",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "622": {
        "name": "Shigella dysenteriae",
        "rank": "species",
        "count_taxon": 0,
        "count_descendant": 338,
        "lineage": {
            "species": "Shigella dysenteriae",
            "genus": "Shigella",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "1224": {
        "name": "Proteobacteria",
        "rank": "phylum",
        "count_taxon": 0,
        "count_descendant": 2875,
        "lineage": {
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "1236": {
        "name": "Gammaproteobacteria",
        "rank": "class",
        "count_taxon": 0,
        "count_descendant": 2875,
        "lineage": {
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "83333": {
        "name": "Escherichia coli K-12",
        "rank": "strain",
        "count_taxon": 0,
        "count_descendant": 385,
        "lineage": {
            "strain": "Escherichia coli K-12",
            "species": "Escherichia coli",
            "genus": "Escherichia",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "91347": {
        "name": "Enterobacterales",
        "rank": "order",
        "count_taxon": 0,
        "count_descendant": 2875,
        "lineage": {
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "131567": {
        "name": "cellular organisms",
        "rank": "no rank",
        "count_taxon": 0,
        "count_descendant": 2875,
        "lineage": {}
    },
    "300267": {
        "name": "Shigella dysenteriae Sd197",
        "rank": "strain",
        "count_taxon": 338,
        "count_descendant": 338,
        "lineage": {
            "strain": "Shigella dysenteriae Sd197",
            "species": "Shigella dysenteriae",
            "genus": "Shigella",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    },
    "511145": {
        "name": "Escherichia coli str. K-12 substr. MG1655",
        "rank": "no rank",
        "count_taxon": 385,
        "count_descendant": 385,
        "lineage": {
            "strain": "Escherichia coli K-12",
            "species": "Escherichia coli",
            "genus": "Escherichia",
            "family": "Enterobacteriaceae",
            "order": "Enterobacterales",
            "class": "Gammaproteobacteria",
            "phylum": "Proteobacteria",
            "superkingdom": "Bacteria"
        }
    }
}
```

## CSV

**Rows**: Taxons

**Columns**:

- `TAXID`: NCBI taxonomy ID
- `name`: Name of the taxon
- `rank`: Taxonomic rank
- `count_taxon`: number of reads mapping to the taxon
- `count_descendant`: number of reads belonging to the descendants of the taxon
- `lineage`: Taxonomic lineage of this taxon, each taxonomic level being separated by a `-` sign.

```
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TAXID  | name                                      | rank         | count_taxon | count_descendant | lineage                                                                                                                                                                                                                           |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1      | root                                      | no rank      | 0           | 2875             |                                                                                                                                                                                                                                   |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 131567 | cellular organisms                        | no rank      | 0           | 2875             |                                                                                                                                                                                                                                   |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2      | Bacteria                                  | superkingdom | 0           | 2875             | superkingdom: Bacteria                                                                                                                                                                                                            |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1224   | Proteobacteria                            | phylum       | 0           | 2875             | phylum: Proteobacteria || superkingdom: Bacteria                                                                                                                                                                                  |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1236   | Gammaproteobacteria                       | class        | 0           | 2875             | class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                                                                                                                                    |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 91347  | Enterobacterales                          | order        | 0           | 2875             | order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                                                                                                         |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 543    | Enterobacteriaceae                        | family       | 2152        | 2875             | family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                                                                           |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 561    | Escherichia                               | genus        | 0           | 385              | genus: Escherichia || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                                                     |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 562    | Escherichia coli                          | species      | 0           | 385              | species: Escherichia coli || genus: Escherichia || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                        |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 83333  | Escherichia coli K-12                     | strain       | 0           | 385              | strain: Escherichia coli K-12 || species: Escherichia coli || genus: Escherichia || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria       |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 511145 | Escherichia coli str. K-12 substr. MG1655 | no rank      | 385         | 385              | strain: Escherichia coli K-12 || species: Escherichia coli || genus: Escherichia || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria       |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 620    | Shigella                                  | genus        | 0           | 338              | genus: Shigella || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                                                        |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 622    | Shigella dysenteriae                      | species      | 0           | 338              | species: Shigella dysenteriae || genus: Shigella || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria                                       |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 300267 | Shigella dysenteriae Sd197                | strain       | 338         | 338              | strain: Shigella dysenteriae Sd197 || species: Shigella dysenteriae || genus: Shigella || family: Enterobacteriaceae || order: Enterobacterales || class: Gammaproteobacteria || phylum: Proteobacteria || superkingdom: Bacteria |
+--------+-------------------------------------------+--------------+-------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
```

## BAM

> Only  generated when running `sam2lca analyze` with the `-b`/`--bam_out` flag

The input alignment file is written as a `bam` file, with the following extra tags:

- `XT` (of type `int`/`i`) set to the TAXID of the LCA assigned to the read
- `XN` (of type `string`/`Z`) set to the scientific name of the LCA assigned to the read
- `XR` (of type `string`/`Z`) set to the taxonomic rank of the LCA assigned to the read

```bam
escherichia_coli_180 355 NC_000913.3 38 1 68M = 148 186 GTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAA DFFAF?DDHAFEBFHGEHIIIGFBFECBFGDBDF?G@HED?FHGHGE>=;@;@@=D@:5:.;;>:@CC AS:i:0 XS:i:0 XM:i:0 XO:i:0 XG:i:0 NM:i:0 MD:Z:68 YS:i:0 YT:Z:CP XT:i:543 XN:Z:Enterobacteriaceae XR:Z:family
shigella_dysenteriae_504 147 NC_007607.1 181065 255 76M = 181033 -108 TGATGACAATTTATTGTCTTATCGTTGTTCTTATGGAACGCTTTTCTGATTGATTTCATATTGGCGAGAGAACAAG @CC>CCCE@EGECHGGGEHEFCIGGGHDFIIIHIIIGJJIIJIJIJJJIJIGEHGEHGJIJJIHF@HGHHHHFDBF AS:i:0 XM:i:0 XO:i:0 XG:i:0 NM:i:0 MD:Z:76 YS:i:0 YT:Z:CP XT:i:300267 XN:Z:Shigella dysenteriae Sd197 XR:Z:strain
```

Reads belonging with these tags can be filtered with [`samtools view`](http://www.htslib.org/doc/samtools-view.html) like this: `samtools view --tag [tag_name]:[value_to_filter] [YOURFILE.bam]`

For Example:

- Reads with the LCA's TAXID equal to `300267`:  `samtools view --tag XT:300267 aligned.sorted.bam`
- Reads with the LCA's rank at `strain` level: `samtools view --tag XR:genus aligned.sorted.sam2lca.bam`
- Reads with the LCA's scientific name being `Shigella dysenteriae Sd197`: `samtools view --tag XN:"Shigella dysenteriae Sd197" aligned.sorted.sam2lca.bam`
