# Tutorial

**Using [sam2lca](https://github.com/maxibor/sam2lca) to identify a plant taxon from `fastq` sequencing files.**

In this tutorial, we'll use the [Angiosperms353](https://academic.oup.com/sysbio/article/68/4/594/5237557) plant markers database to identify a plant species present in our sequencing data. The Angiosperms353 database consists of up to 353 universal Angiosperms (flowering plants) gene markers that are derived from the [1000 plant transcriptomes](https://www.nature.com/articles/s41586-019-1693-2) project.

## Installing all tools for this tutorial

For this tutorial, a dedicated conda-environment is available to ease the reproducibility.

Download the environment:

```bash
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/environment.yaml
```

Installing and activating environment:

```bash
conda env create -f environment.yaml
conda activate sam2lca_tutorial
```

## Getting the reference database

First we need to download and decompress the reference `fasta` database, which in this case, has already been pre-formatted for sam2lca.

> For the sake of this tutorial, we use a reduced version of the [angiosperms353](https://academic.oup.com/sysbio/article/68/4/594/5237557) marker set. The full version is available for download, see [sam2lca.readthedocs.io/en/latest/databases.html](https://sam2lca.readthedocs.io/en/latest/databases.html).  

```bash
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/tutorial_db.fa.gz
gunzip tutorial_db.fa.gz
```

## Indexing the database with Bowtie2

In this tutorial, we're going to work with the read aligner [Bowtie2](https://github.com/BenLangmead/bowtie2), but other aligners like [BWA](http://bio-bwa.sourceforge.net/) work also just fine.

Before being able to do any alignment, we need to index the Angiosperms353 database with Bowtie2:

```bash
bowtie2-build tutorial_db.fa angiosperms353
```

## Preparing `fastq` sequencing files

> This step might be a bit long, especially if you have many references present in your database. You may want to speed it up by parallelizing it using the `--threads` option.

Prior to aligning the sequencing data, we need to download and process the sequencing data, e.g. to remove adapter sequences. 

Downloading the paired-end DNA sequencing compressed `fastq` files

```bash
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/metagenome.1.fastq.gz
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/metagenome.2.fastq.gz
```

- Performing adapter-clipping and quality trimming with [fastp](https://github.com/OpenGene/fastp)

```bash
$ fastp -i metagenome.1.fastq.gz -I metagenome.2.fastq.gz -o metagenome_trimmed.R1.fastq.gz -O metagenome_trimmed.R2.fastq.gz
```

## Alignment with Bowtie2

After having prepared both the Angiosperms353 database and the FastQ sequencing files, we now align the sequencing data against the references using BowTie2.

> The important aspect here is to allow that multiple alignments can be reported for each read to ensure that all potential hits are reported. This is done by using the `-a` flag of Bowtie2 for reporting alignments, or `-k` for reporting up to *N* alignments.

Here, we will allow the reporting of up to 50 alignments per read.

```bash
bowtie2 -x angiosperms353 -k 50 -1 metagenome_trimmed.R1.fastq.gz -2 metagenome_trimmed.R2.fastq.gz | samtools sort -O bam > metagenome.sorted.bam
samtools index metagenome.sorted.bam
```

## Optional but (highly) recommended: [bamAlignCleaner](https://github.com/maxibor/bamAlignCleaner)

Like many other tools working with `bam`/`cram` files, [sam2lca](https://github.com/maxibor/sam2lca) relies on the file index to facilitate a rapid and parallel access to the aligned segments.

However, because we aligned our sequencing data against a database containing (most likely) a lot more reference sequences than what we can potentially align our reads to, the index of the file, based on the `bam`/`cram` header is huge and will slow down the [I/O](https://en.wikipedia.org/wiki/Input/output) by many folds.

To circumvent this, we advise you to run [bamAlignCleaner](https://github.com/maxibor/bamAlignCleaner) on your alignment files, and then reindexing them, before processing them further with sam2lca.

> For the sake of the tutorial, we use a smaller `fasta` reference sequence database. This step is therefore not really necessary in this tutorial.

```bash
bamAlignCleaner metagenome.sorted.bam | samtools sort > metagenome.cleaned.sorted.bam
samtools index metagenome.cleaned.sorted.bam
```

## Running sam2lca

Once we have our alignment file, here in `bam` format, we can now run [sam2lca](https://github.com/maxibor/sam2lca) to identify which plants shed some of its DNA in our sequencing file.

First, we need to set up the sam2lca acc2tax database for *plant markers*, with NCBI taxonomic identifiers. This step downloads and creates the taxonomy and accession to taxid conversion databases. These databases allow sam2lca to convert the accession ids of the reference genome sequences to their respective taxonomic ids and prepares this information to be accessible within sam2lca.

```bash
sam2lca update-db --taxonomy ncbi --acc2tax plant_markers
```

Let's check which sam2lca databases are now available:

```bash
sam2lca list-db
```

Finally, we run sam2lca with the *plant markers* database.

> To make sure that we don't accidentally run the LCA algorithm on DNA sequences that are unlikely to belong to the same clade, we will only run the LCA for all references aligned to each read that have a identity greater than 90%. Depending on the type of database, you might want to adjust the sequence identity threshold or try different ones.

```bash
sam2lca analyze --acc2tax plant_markers -b -i 0.9 metagenome.cleaned.sorted.bam
```

Let's look at the results that are summarized in the file `metagenome.sorted.sam2lca.csv`

We see that the only species present in our data, is [*Cannabis sativa*](https://en.wikipedia.org/wiki/Cannabis_sativa), which is indeed what was present in our sample ! (it was a simulated dataset 😉)

```
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TAXID   | name               | rank         | count_taxon | count_descendant | lineage                                                                                                                                                                                                                                      |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3398    | Magnoliopsida      | class        | 5           | 79               | class: Magnoliopsida || clade: Embryophyta || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 58024   | Spermatophyta      | clade        | 0           | 79               | clade: Embryophyta || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                 |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 131567  | cellular organisms | no rank      | 0           | 79               |                                                                                                                                                                                                                                              |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2759    | Eukaryota          | superkingdom | 0           | 79               | superkingdom: Eukaryota                                                                                                                                                                                                                      |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 33090   | Viridiplantae      | kingdom      | 0           | 79               | kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                                                                                            |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 35493   | Streptophyta       | phylum       | 0           | 79               | phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                                                                    |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 131221  | Streptophytina     | subphylum    | 0           | 79               | subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                                       |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3193    | Embryophyta        | clade        | 0           | 79               | clade: Embryophyta || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                 |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 58023   | Tracheophyta       | clade        | 0           | 79               | clade: Embryophyta || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                 |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 78536   | Euphyllophyta      | clade        | 0           | 79               | clade: Embryophyta || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                                                 |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1       | root               | no rank      | 0           | 79               |                                                                                                                                                                                                                                              |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1437183 | Mesangiospermae    | clade        | 0           | 74               | clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 71240   | eudicotyledons     | clade        | 0           | 74               | clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 91827   | Gunneridae         | clade        | 0           | 74               | clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1437201 | Pentapetalae       | clade        | 2           | 74               | clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 71275   | rosids             | clade        | 1           | 72               | clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 91835   | fabids             | clade        | 3           | 71               | clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                                         |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3744    | Rosales            | order        | 6           | 68               | order: Rosales || clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                                       |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3481    | Cannabaceae        | family       | 21          | 62               | family: Cannabaceae || order: Rosales || clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                                                |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3482    | Cannabis           | genus        | 0           | 41               | genus: Cannabis || family: Cannabaceae || order: Rosales || clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota                             |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3483    | Cannabis sativa    | species      | 41          | 41               | species: Cannabis sativa || genus: Cannabis || family: Cannabaceae || order: Rosales || clade: Embryophyta || class: Magnoliopsida || subphylum: Streptophytina || phylum: Streptophyta || kingdom: Viridiplantae || superkingdom: Eukaryota |
+---------+--------------------+--------------+-------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
```

> Note that out of all the reads in our sample, with the sam2lca parameters we used, 87.8% were classified (79 out of 90 trimmed reads). However, only 41 of them were assigned at the species level. This is due to the nature of these angiospersm353 markers: they are gene markers, hence relatively highly conserved among plants, which explains why the LCA is bringing many reads to a lower resolution taxonomic level.

Last but not least, we can also have a look at the `bam` file created by sam2lca: `metagenome.cleaned.sorted.sam2lca.bam`

```bash
$ samtools view metagenome.cleaned.sorted.sam2lca.bam | head -n 2
oneKP_189 329 5427_Ryparosa_wrayi 621 255 76M = 621 0 GGAGATTCAGAAGATGGCAAAATCAATTAAGGAACTGAAGAAGGAAAATTCATTCTTGAAGAGCAAGACTGAGAAA D=BFCFCFHGIFCDGHBFH<@H?CFHIGHGIIGHDGEFGF@<=FHICHEIGHHHHFFFCECB>?;;ACDCCCDCDC AS:i:-46 XM:i:9 XO:i:0XG:i:0 NM:i:9 MD:Z:0A2A3G28C1G10A11G5T0G7 YT:Z:UP XT:i:12908 XN:Z:unclassified sequences
oneKP_95 329 5921_Pappea_capensis 2336 255 76M = 2336 0 AGATTCAGTCTGATATTGACAAAAATAGTACCGACATCAACCGTCACAAGGTTCAAATAGAGACTGCCCAGAAAAT DBDFHHFHIIIGIBHFHIIIIIIIGIIIIIIIIIIIIIIIIAEHIIFGDHHHFFDBCECCCB@BCCCCCCCBCCBC AS:i:-43 XM:i:8 XO:i:0XG:i:0 NM:i:8 MD:Z:28C2T2G14C14G1G0A2A5 YT:Z:UP XT:i:12908 XN:Z:unclassified sequences
oneKP_76 329 6270_Parartocarpus_venenosus 169 255 76M = 169 0 GACATAGACTATGGGAATGATGTGTTAACCTTGAAGCTTGGTGATTTAGGAACGTATGTGTTGAACAAACAAACTC ?D;DBCBDAE:AEEAED?<C?<D?9?D?DDD?DADCC=.-;'>?77?).6;6>>A5-;;>>:;<5=1>;:2)(89> AS:i:-35 XM:i:9XO:i:0 XG:i:0 NM:i:9 MD:Z:26G2G0G10G5G5A5T0A1A13 YT:Z:UP XT:i:12908 XN:Z:unclassified sequences
oneKP_189 329 5427_Seidelia_firmula 609 255 76M = 609 0 GGAGATTCAGAAGATGGCAAAATCAATTAAGGAACTGAAGAAGGAAAATTCATTCTTGAAGAGCAAGACTGAGAAA D=BFCFCFHGIFCDGHBFH<@H?CFHIGHGIIGHDGEFGF@<=FHICHEIGHHHHFFFCECB>?;;ACDCCCDCDC AS:i:-46 XM:i:9 XO:i:0XG:i:0 NM:i:9 MD:Z:0A6G10G2G5C8C12A0T0G24 YT:Z:UP XT:i:12908 XN:Z:unclassified sequences
oneKP_95 329 5921_Placodiscus_turbinatus 980 255 76M = 980 0 AGATTCAGTCTGATATTGACAAAAATAGTACCGACATCAACCGTCACAAGGTTCAAATAGAGACTGCCCAGAAAAT DBDFHHFHIIIGIBHFHIIIIIIIGIIIIIIIIIIIIIIIIAEHIIFGDHHHFFDBCECCCB@BCCCCCCCBCCBC AS:i:-43 XM:i:8XO:i:0 XG:i:0 NM:i:8 MD:Z:28C2T2G14A14A1G0A2A5 YT:Z:UP XT:i:12908 XN:Z:unclassified sequences
oneKP_144 73 6406_Cannabis_sativa 103 255 76M = 103 0 TCGAATACTGGACGGCTTTTCTTCGAACGAGCTATTGGAGCTTTTAGTATTGCCGAAATGCTGTTTACCTCTAGTC FFFGHHHHJJJJJJJJJJJJJHIJJIJJJJJIJJJJGIJJJJJJJJJJIIJJJJIIJIHFGGGH;DHCAEHFFDFB AS:i:0 XM:i:0 XO:i:0 XG:i:0NM:i:0 MD:Z:76 YT:Z:UP XT:i:3483 XN:Z:Cannabis sativa XR:Z:species
```

Note the `XT`, `XN`, and `XR` tags at the end of each line. The value of these tags (the value  after the `:` ) is set by sam2lca to be respectively, the taxonomy IDs, the taxon name, and the taxon rank attributed to each sequencing read by the LCA algorithm.

To filter for these tags, we can make use of the filter expressions included in `samtools view`. For example, to obtain all reads assigned to the taxonomic level *class*, we can use `samtools view -e '[XR] == "class"' metagenome.cleaned.sorted.sam2lca.bam` or `samtools view --tag XR:class metagenome.cleaned.sorted.sam2lca.bam`.

Also note that for the first 5 reads (`oneKP_189`, `oneKP_95`, `oneKP_76`, `oneKP_189`), sam2lca classified them as "unclassified sequences" because they didn't pass the sequence identity threshold of 90%. Only the reads `oneKP_144` aligned to `6406_Cannabis_sativa` has an identity high enough to be classified by sam2lca.
