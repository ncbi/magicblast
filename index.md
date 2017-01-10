---
layout: page
title: "NCBI Magic-BLAST Documentation"
---

Magic-BLAST is a tool for mapping large next-generation RNA or DNA sequencing
runs against a whole genome or transcriptome. Each alignment optimizes
a composite score, taking into account simultaneously the two reads of
a pair, and in case of RNA-seq, locating the candidate introns and adding
up the score of all exons. This is very different from other versions of
BLAST, where each exon is scored as a separate hit and read-pairing is
ignored.

[Download Magic-BLAST from here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST)

Magic-BLAST incorporates within the NCBI BLAST code framework ideas
developed in the NCBI Magic pipeline, in particular hit extensions by
local walk and jump [(http://www.ncbi.nlm.nih.gov/pubmed/26109056)](http://www.ncbi.nlm.nih.gov/pubmed/26109056), and recursive clipping of
mismatches near the edges of the reads, which avoids accumulating
artefactual mismatches near splice sites and is needed to distinguish
short indels from substitutions near the edges.

We call the whole next generation run (from Illumina, Roche-454, ABI, or
another sequencing platform excluding SOLiD), a query. The input reads may
be provided as SRA accession or a file in a SRA, FASTA, and FASTQ format.
Read pairs can be presented as parallel files, or as successive reads in a
single file.

The reference genome or transcriptome can be given as a BLAST database
or a FASTA file. It is preferable to use BLAST database for large genomes,
such as human, or transcript collections, such as all of RefSeq, Ensembl,
or AceView. See here on [how to create a BLAST database](/magicblast/doc/blastdb.html).

The full list of options is listed when you use -help option.

Thank you for testing this code and providing us with feedback. Please,
let us know of any desired option, problem or difficulty.

E-mail boratyng@ncbi.nlm.nih.gov or blast-help@ncbi.nlm.nih.gov with questions or comments.
