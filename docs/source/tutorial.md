# Tutorial

**Using [sam2lca](https://github.com/maxibor/sam2lca) to identify a plant taxon from `fastq` sequencing files.**

In this tutorial, we'll use the [Angiosperms353](https://academic.oup.com/sysbio/article/68/4/594/5237557) plant markers database to identify a plant species present in our sequencing data. The Angiosperms353 database consists of up to 353 universal Angiosperms (flowering plants) gene markers that are derived from the [1000 plant transcriptomes](https://www.nature.com/articles/s41586-019-1693-2) project.

## Installing all tools for this tutorial

For this tutorial, a dedicated conda-environment is available to ease the reproducibility.

Download the environment:

```bash
$ wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/environment.yaml
```

Installing and activating environment:

```bash
$ conda env create -f environment.yaml
$ conda activate sam2lca_tutorial
```


## Getting the reference database

First we need to download and decompress the Angiosperms353 database, which has already been pre-formatted for sam2lca

```bash
$ wget -O angiosperms353_markers.fa.gz https://edmond.mpdl.mpg.de/api/access/datafile/101862
$ gunzip angiosperms353_markers.fa.gz
```

## Indexing the database with Bowtie2

In this tutorial, we're going to work with the read aligner [Bowtie2](https://github.com/BenLangmead/bowtie2), but other aligners like [BWA](http://bio-bwa.sourceforge.net/) work also just fine.

Before being able to do any alignment, we need to index the Angiosperms353 database with Bowtie2:

```bash
$ bowtie2-build angiosperms353_markers.fa angiosperms353
```

## Preparing the `fastq` sequencing files

Prior to aligning the sequencing data, we need to download and process the sequencing data, e.g. to remove adapter sequences. 

Downloading the paired-end DNA sequencing compressed `fastq` files

```bash
$ wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/metagenome.1.fastq.gz
$ wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/metagenome.2.fastq.gz
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
$ bowtie2 -x angiosperms353 -k 50 -1 metagenome_trimmed.R1.fastq.gz -2 metagenome_trimmed.R2.fastq.gz | samtools sort -O bam > metagenome.sorted.bam
$ samtools index metagenome.sorted.bam
```

## Running sam2lca

Once we have our alignment file, here in `bam` format, we can now run [sam2lca](https://github.com/maxibor/sam2lca) to identify which plants shed some of its DNA in our sequencing file.

First, we need to set up the sam2lca database for *plant markers*. This steps matches the accession ids of the reference genomes with their respective taxonomic ids and prepares this information to be accessible within sam2lca.

```bash
$ sam2lca -m plant_markers update-db
```

Finally, we run sam2lca with the *plant markers* database. 

> To make sure that we don't accidentally run the LCA algorithm on DNA sequences that are unlikely to belong to the same clade, we will only run the LCA for all references aligned to each read that have a identity greater than 90%. Depending on the type of database, you might want to adjust the sequence identity threshold or try different ones.

```bash
$ sam2lca -m plant_markers analyze -i 0.9 metagenome.sorted.bam
```

Let's look at the results, for example the file `metagenome.sorted.sam2lca.csv`

We see that the only species present in our data, is [*Cannabis sativa*](https://en.wikipedia.org/wiki/Cannabis_sativa), which is indeed what was present in our sample ! (it was a simulated dataset 😉)



```
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TAXID | name            | rank    | count | lineage                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3483  | Cannabis sativa | species | 38    | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids' - 'order': 'Rosales' - 'family': 'Cannabaceae' - 'genus': 'Cannabis' - 'species': 'Cannabis sativa' |
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3481  | Cannabaceae     | family  | 34    | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids' - 'order': 'Rosales' - 'family': 'Cannabaceae'                                                      |
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3744  | Rosales         | order   | 9     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids' - 'order': 'Rosales'                                                                                |
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3398  | Magnoliopsida   | class   | 5     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida'                                                                                                                                                                                                                                                        |
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 91835 | fabids          | clade   | 2     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids'                                                                                                     |
+-------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
```


> Note that out of all the reads in our sample, 98.9% were classified. However, only 38 of them (out of 88) were assigned at the species level. This is due to the nature of these angiospersm353 markers: they are gene markers, hence relatively highly conserved among plants, which explains why the LCA is bringing many reads to a lower resolution taxonomic level.
