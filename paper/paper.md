---
title: 'sam2lca: Lowest Common Ancestor for SAM/BAM/CRAM alignment files'
tags:
  - Python
  - metagenomics
  - microbiome
  - alignment
  - sam
  - bam
  - cram
  - htslib
  - taxonomy
  - LCA
  - ancestor
authors:
  - name: Maxime Borry # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0001-9140-7559
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Alexander Hübner # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0003-3572-9996
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Christina Warinner
    orcid: 0000-0002-4528-5877
    affiliation: "1, 2, 3"
affiliations:
 - name: Microbiome Sciences Group, Max Planck Institute for Evolutionary Anthropology, Department of Archaeogenetics, Leipzig, Germany
   index: 1
 - name: Faculty of Biological Sciences, Friedrich-Schiller Universität Jena, Jena, Germany
   index: 2
 - name: Department of Anthropology, Harvard University, Cambridge, MA, United States of America
   index: 3
date: 4 April 2022
bibliography: paper.bib

---

# Summary

sam2lca is a program performing reference sequence disambiguation for reads mapping to multiple reference sequences in a shotgun metagenomics sequencing dataset. To do so, it takes as input the common SAM sequence alignment format and applies the lowest common ancestor algorithm.

# Statement of need

The rapidly decreasing cost of massively parallel short-read DNA sequencing technologies has enabled the genetic characterization of entire ecological communities, a technique known as shotgun metagenomics.

In a typical shotgun metagenomics approach, after the DNA of an ecological community has been sequenced, it is compared to a genetic reference database of organisms with known taxonomy. Even though the number of DNA sequences and genomes in reference databases is constantly growing, there are still instances where a query sequence will not have a direct match in a reference database, and it will instead weakly align to one or more distantly related reference organisms. Furthermore, when analyzing short DNA sequences, a query DNA sequence will often match equally well to more than one reference organism, posing a challenge for its taxonomic assignation.

One solution to this problem is to apply a lowest common ancestor algorithm (LCA) (\autoref{fig:fig1}) during taxonomic profiling to place such ambiguous assignments higher in a taxonomic tree, where they can be more confidently assigned. This idea was first implemented for metagenomics with the MEGAN program [@huson2007megan].

Many programs have since been developed to perform LCA during taxonomic profiling. For example, MALT [@herbig2017malt] and MetaPhlan [@segata2012metagenomic] perform LCA and taxonomic profiling after DNA sequence alignment, while other programs, such as Kraken2 [@kraken2] and Centrifuge [@kim2016centrifuge], are alignment-free methods that apply LCA after k-mer matching. While combining the steps of database matching and LCA into one program can be useful, it also limits user choice for the selection of different alignment or k-mer matching programs.

With sam2lca, we propose to decouple the LCA step from the alignment step to allow the end-user to freely choose from one of the many DNA sequence aligner programs available, such as Bowtie2 [@bowtie2], bwa [@bwa], bbmap [@bbmap], or minimap2 [@minimap2]. Each of these aligners exports the sequence alignments in the widely adopted Sequence Alignment Map format (SAM) [@sam_format], or in its binary (BAM), or compressed representation (CRAM), which sam2lca uses as an input.

The use of the SAM file format enables easier integration of sam2lca in a wide variety of analysis workflows, which often already contain steps generating or using SAM/BAM/CRAM files, and allows for an easy subsequent analysis using well-established programs, such as SAMtools [@sam_format].

![Example of the LCA algorithm with NCBI TAXIDs. Taxons and their LCA are displayed in the same color. The lineage for each taxon is shown with a one letter code for the rank, and the corresponding TAXID. The LCA of s:562 (*E. coli* species) and s:622 (*S. dysenteriae* species) is f:543 (*Enterobacteriaceae* family). The LCA of s:82981 (*L. grimontii* species) and s:158841 (*L. richardii* species) is g:82980 (*Leminorella* genus). The LCA of s:2562891 (*E. alba* species), s:623 (*S. flexneri* species) and s:2498113 (*J. zhutongyuii* species) is o:91347 (Enterobacterales order) \label{fig:fig1}](figures/figure1.png)

# Implementation

sam2lca is a program written in Python, which takes as an input an indexed and sorted SAM/BAM/CRAM alignment file. Broadly, the program consists of four main steps. First, reference sequence accessions, present in the BAM file header section, are converted to taxonomic identifiers (TAXID) using a RocksDB persistent key-value store [@dong2021rocksdb]. The alignment section of the BAM file is then parsed with Pysam [@pysam] and a dictionary is created to match single and multi-mapping query sequences/reads to the TAXID(s) of their matching  reference sequence(s). Next, if a read has been matched to multiple TAXIDs, the LCA implementation of Taxopy [@taxopy] is used to attribute it to the lowest common ancestor, using the NCBI taxonomy by default. Finally, each TAXID is used to retrieve its associated taxon's scientific name and taxonomic lineage, and results are saved in a JSON and CSV file. Optionally, a BAM file, similar to the input file, can be generated. This BAM file contains for each read an additional XT tag added to report the TAXID of the LCA for each read, an XN tag for the taxon's scientific name, and finally an XR tag for the taxon's rank. sam2lca is distributed through pip and conda, and the documentation and tutorials are available at [sam2lca.readthedocs.io](https://sam2lca.readthedocs.io)

# Acknowledgements

This research was supported by the Werner Siemens Stiftung (M.B. and C.W.) and by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy EXC 2051, Project-ID 390713860 (A.H. and C.W.).

# References
