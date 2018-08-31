#!/bin/tcsh
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
#  Author: Jean Thierry-Mieg, Greg Boratyn
#

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


