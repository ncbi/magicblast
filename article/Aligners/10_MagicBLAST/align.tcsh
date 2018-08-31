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
#  Author: Jean Thierry-Mieg
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

# create a BLAST database if not present
if (! -e $genomeDir/genome.nhr) then
  echo $genomeDir/genome.nhr
  gunzip -c Fasta/$run/genome.gz > $genomeDir/my_genome.fasta
  bin/makeblastdb -in $genomeDir/my_genome.fasta -out $genomeDir/genome -dbtype nucl -parse_seqids
  \rm $genomeDir/my_genome.fasta
endif

mkdir -p $mm/$run

# run magicblast
echo "$mm $run $reads $mates"
ls -ls  $genomeDir/genome.nhr
echo "$reads"
ls -ls $reads
echo "$reads"
set infmt=`echo $reads | gawk '/fastq/{print "fastq";next}/fasta/{print "fasta"}'`


if (-e  $genomeDir/genome.nhr && -e $reads && ! -e $mm/$run/$mm.$run.sam) then
  set mmm="-query_mate $mates"
  if (X$mates == X) set mmm=""
  echo "time bin/magicblast -query $reads $mmm  -infmt $infmt -db $genomeDir/genome -num_threads $nThreads "
  time bin/magicblast -query $reads $mmm -infmt $infmt -db $genomeDir/genome -num_threads $nThreads > $out.sam_unsorted
  ls -ls $out.sam_unsorted
  time sort $out.sam_unsorted >  $out.sam_sorted
  gzip $out.sam_sorted
  ls -ls $out.sam_sorted.gz
  # rm   $mm/$run/$mm.$run.sam_unsorted
else
  echo "Did not find $genomeDir/genome.nhr or  $reads or found $mm/$run/$mm.$run.sam"
  ls -ls  $mm/$run/$mm.$run.sam
endif


