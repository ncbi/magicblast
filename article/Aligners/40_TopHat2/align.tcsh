#!/bin/tcsh

set mm=$1
set run=$2
set genome=$3
set reads=$4
set mates=$5

set genomeF=Fasta/$run/genome.gz
set out=$mm/$run/$mm.$run

set nThreads=4
set genomeDir=Aligners/$mm/$genome
if (! -e $genomeDir) mkdir -p $genomeDir

mkdir -p $genomeDir
set path=($path `pwd`/bin/bowtie2)

# create an index if not present
if (! -e $genomeDir/genome.1.bt2) then
  echo $genomeDir/genome.1.bt2
  gunzip -c Fasta/$run/genome.gz > $genomeDir/my_genome.fasta
  bin/bowtie2/bowtie2-build $genomeDir/my_genome.fasta $genomeDir/genome
  \rm $genomeDir/my_genome.fasta
endif

mkdir -p $mm/$run

# run tophat
echo "$mm $run $reads $mates"
ls -ls  $genomeDir/genome.1.bt2
echo "$reads"
ls -ls $reads
echo "$reads"

if (-e  $genomeDir/genome.1.bt2 && -e $reads && ! -e $mm/$run/$mm.$run.sam) then
  echo "time bin/tophat2/tophat2 -p $nThreads -o ${out}_dir $genomeDir/genome $reads $mates "
  time bin/tophat2/tophat2 -p $nThreads -o ${out}_dir $genomeDir/genome $reads $mates
  samtools view -h ${out}_dir/accepted_hits.bam >$out.sam_unsorted 
  time sort  $mm/$run/$mm.$run.sam_unsorted >  $mm/$run/$mm.$run.sam_sorted
  gzip  $mm/$run/$mm.$run.sam_sorted
  \rm   $mm/$run/$mm.$run.sam_unsorted
else
  echo "Did not find $genomeDir/genome.1.ht2 or  $reads or found $mm/$run/$mm.$run.sam"
  ls -ls  $mm/$run/$mm.$run.sam
endif


