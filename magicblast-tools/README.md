Magic-BLAST tools
===

A a few scripts to scripts to postprocess SAM/BAM files.


Features
--------
* get-introns.py - collect intron locations.
* get-transcripts.py - assemble transcript sequences from a genome using a GFF or GTF annotation.
* combine-genome-transcript.py - iterate over read alignments to a genome and transcripts, select better scoring alignments, remap transcript alignments to the genome, and save them in a SAM or BAM file.


## Dependencies
The programs in this directory work with Python 3.6+ and require these packages:
* pysam
* pyfaidx
* pandas

To install them run:
```
pip install -r requirements.txt
```


Usage
-----

### Get intron locations from a SAM/BAM file
```
get-introns.py --bam <SAM/BAM file> --gff <GFF/GTF/SAM/BAM file with annotation> --introns <output file>
```

The output is a tab delimited file with intron locations marked as KNOWN or NEW.


### Assemble transcript sequences from a genome
```
get-transcripts.py --genome <FASTA genome> --gff <GFF/GTF file>
```
The transcripts will be written to the standard output.


### Combine alignments to genome and transcripts
```
combine-genome-transcripts.py --to-genome <SAM/BAM genome alignments> --to-transcripts <SAM/BAM transcript alignments> --gff <GFF/GTF annotation> --out <output SAM/BAM file> [-b]
```
