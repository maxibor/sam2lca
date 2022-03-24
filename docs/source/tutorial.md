# Tutorial

**Using [sam2lca](https://github.com/maxibor/sam2lca) to identify a plant taxon from a `fastq` sequencing files.**

In this tutorial, we'll use the [Angiosperms353](https://academic.oup.com/sysbio/article/68/4/594/5237557) plant markers database, which consists of up to 353 universal Angiosperms (flowering plants) gene markers, derived from the [1000 plant transcriptomes](https://www.nature.com/articles/s41586-019-1693-2) project, to identify a the plant species present in our sequencing data.

## Installing all tools for this tutorial

For this tutorial, a dedicated conda-environment is available to ease the reproducibility.

Download environment

```bash
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/environment.yaml
```

Installing and activating environment

```bash
conda env create -f environment.yaml
conda activate sam2lca_tutorial
```

## Getting the reference database

First we need to download and decompress the Angiosperms353 database, which has already been pre-formatted for sam2lca

```bash
wget -O angiosperms353_markers.fa.gz https://edmond.mpdl.mpg.de/api/access/datafile/101862
gunzip angiosperms353_markers.fa.gz
```

## Indexing the database with Bowtie2

In this tutorial, we're going to work with the [Bowtie2](https://github.com/BenLangmead/bowtie2) read aligner, but other aligners like [BWA](http://bio-bwa.sourceforge.net/) are also just fine.

Before doing any alignment, we need to index the angiosperms353 database with bowtie2

```bash
bowtie2-build angiosperms353_markers.fa angiosperms353
```

> This step might be a bit long, because of the many references present in this database, you may want to speed it up by parallelizing it using the `--threads` option.

## Preparing `fastq` sequencing files

- Downloading the paired-end DNA sequencing compressed `fastq` files

```bash
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/metagenome.1.fastq.gz
wget https://raw.githubusercontent.com/maxibor/sam2lca/master/docs/tutorial/data/metagenome.2.fastq.gz
```

- Performing adapter-clipping and quality trimming with [fastp](https://github.com/OpenGene/fastp)

```bash
fastp -i metagenome.1.fastq.gz -I metagenome.1.fastq.gz -o metagenome_trimmed.R1.fastq.gz -O metagenome_trimmed.R2.fastq.gz
```

## Alignment with Bowtie2

> The important aspect here is to allow multiple alignments reporting for each read, to make sure that all potential hits are reported. This is done by using the `-a` flag of bowtie for reporting alignemnts, or `-k` for reporting up to *N* alignments.

Here, we will allow the reporting of up to 50 alignments per read.

```bash
bowtie2 -x angiosperms353 -k 50 -1 metagenome_trimmed.R1.fastq.gz -2 metagenome_trimmed.R2.fastq.gz | samtools sort -O bam > metagenome.sorted.bam
samtools index metagenome.sorted.bam
```

## Optional but (highly) recommended: [bamAlignCleaner](https://github.com/maxibor/bamAlignCleaner)

Like many other tools working with `bam`/`cram` files, [sam2lca](https://github.com/maxibor/sam2lca) relies on the file index to facilitate a rapid and parallel access to the aligned segments.

However, because we aligned our sequencing data against a database containing (most likely) a lot more reference sequences than what we can potentially align our reads to, the index of the file, based on the `bam`/`cram` header is huge and will slow down the [I/O](https://en.wikipedia.org/wiki/Input/output) by many folds.

To circumvent this, we advise you to run [bamAlignCleaner](https://github.com/maxibor/bamAlignCleaner) on your alignment files, and then reindexing them, before processing them further with sam2lca.

```bash
bamAlignCleaner metagenome.sorted.bam | samtools sort > metagenome.cleaned.sorted.bam
samtools index metagenome.cleaned.sorted.bam
```

## Running sam2lca

Once we have our alignment file, here in `bam` format, we can now run [sam2lca](https://github.com/maxibor/sam2lca) to identify which plants shed some of its DNA in our sequencing file.

First, we need to set up the sam2lca database for *plant markers*

```bash
sam2lca -m plant_markers update-db
```

Finally, run sam2lca with the *plant markers* database.

> To make sure that we don't accidentally run the LCA algorithm on DNA sequences that are unlikely to belong to the same clade, we will only run the LCA for all references aligned to each read that have a identity greater than 90%. Depending on the type of database, you might want to play with sequence identity threshold.

```bash
sam2lca -m plant_markers analyze -b -i 0.9 metagenome.cleaned.sorted.bam
```

Let's look at the results, for example the file `metagenome.sorted.sam2lca.csv`

We see that the only species present in our data, is [*Cannabis sativa*](https://en.wikipedia.org/wiki/Cannabis_sativa), which is indeed what was present in our sample ! (it was a simulated dataset 😉)

```
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TAXID   | name            | rank    | count | lineage                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3483    | Cannabis sativa | species | 41    | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids' - 'order': 'Rosales' - 'family': 'Cannabaceae' - 'genus': 'Cannabis' - 'species': 'Cannabis sativa' |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3481    | Cannabaceae     | family  | 21    | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids' - 'order': 'Rosales' - 'family': 'Cannabaceae'                                                      |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3744    | Rosales         | order   | 6     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids' - 'order': 'Rosales'                                                                                |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3398    | Magnoliopsida   | class   | 5     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida'                                                                                                                                                                                                                                                        |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 91835   | fabids          | clade   | 3     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids' - 'clade': 'fabids'                                                                                                     |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1437201 | Pentapetalae    | clade   | 2     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae'                                                                                                                                             |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 71275   | rosids          | clade   | 1     | 'no rank': 'root' - 'no rank': 'cellular organisms' - 'superkingdom': 'Eukaryota' - 'kingdom': 'Viridiplantae' - 'phylum': 'Streptophyta' - 'subphylum': 'Streptophytina' - 'clade': 'Embryophyta' - 'clade': 'Tracheophyta' - 'clade': 'Euphyllophyta' - 'clade': 'Spermatophyta' - 'class': 'Magnoliopsida' - 'clade': 'Mesangiospermae' - 'clade': 'eudicotyledons' - 'clade': 'Gunneridae' - 'clade': 'Pentapetalae' - 'clade': 'rosids'                                                                                                                         |
+---------+-----------------+---------+-------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
```

> Note that out of all the reads in our sample, with the sam2lca parameters we used, 87.8% were classified (79 out of 90 trimmed reads). However, only 38 of them were assigned at the species level. This is due to the nature of these angiospersm353 markers: they are gene markers, hence relatively highly conserved among plants, which explains why the LCA is bringing many reads to a lower resolution taxonomic level.

Last but not least, we can also have a look at the `bam` file created by sam2lca: `metagenome.cleaned.sorted.sam2lca.bam`

```bash
$ samtools view metagenome.cleaned.sorted.sam2lca.bam | head -n 2
oneKP_189 73 5427_Cannabis_sativa 1092 34 76M = 1092 0 GGAGATTCAGAAGATGGCAAAATCAATTAAGGAACTGAAGAAGGAAAATTCATTCTTGAAGAGCAAGACTGAGAAA D=BFCFCFHGIFCDGHBFH<@H?CFHIGHGIIGHDGEFGF@<=FHICHEIGHHHHFFFCECB>?;;ACDCCCDCDC AS:i:0 XN:i:0 XM:i:0 XO:i:0 XG:i:0 NM:i:0 MD:Z:76 YT:Z:UP XT:i:3398
oneKP_189 133 5427_Cannabis_sativa 1092 0 * = 1092 0 GGAGATTCAGAAGATGGCAAAATCAATTAAGGAACTGAAGAAGGAAAATTCATTCTTGAAGAGCAAGACTGAGAAA D=BFCFCFHGIFCDGHBFH<@H?CFHIGHGIIGHDGEFGF@<=FHICHEIGHHHHFFFCECB>?;;ACDCCCDCDC YT:Z:UP XT:i:3398
```

Note the `XT` tag at the end of each line. The value of this tag (the integer number after `XT:i:`) is set by sam2lca to the NCBI Taxonomy IDs attributed to each sequencing read by the LCA algorithm.
