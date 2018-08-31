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
set target=$3
set reads=$4
set mates=$5

set genome=Fasta/$run/genome.gz
set out=$mm/$run/$mm.$run
mkdir -p $mm/$run
set mm2=21_HISAT2

set nThreads=4
set genomeDir=Aligners/$mm2/$target
if (! -e $genomeDir) mkdir -p $genomeDir

if (! $?TMPDIR) then
  set TMPDIR=/tmp
  if (-d /export/home/TMP) set TMPDIR=/export/home/TMP
endif

# construct the hisat genome index
if (-e $genome  && ! -e $genomeDir/genome.1.ht2) then
  echo "gunzip -c $genome > $genomeDir/mygenome.fasta"
        gunzip -c $genome > $genomeDir/mygenome.fasta
  ls -ls $genomeDir/mygenome.fasta
  if (-e $genomeDir/mygenome.fasta) then
    time bin/hisat2-master/hisat2-build $genomeDir/mygenome.fasta $genomeDir/genome
    touch $genomeDir/done
    # \rm  $genomeDir/mygenome.fasta
  endif
endif

if (-e $genomeDir/genome.1.ht2) then
  echo "HISAT2 index is ready"
else
  echo "missing HISAT2 index $genomeDir/genome.1.ht2"
  goto done
endif

if ($reads == "") then
  echo "missing parametets 4 which should be the reads file"
  goto done
endif

if (! -e $reads) then
  echo "cannot find the reads file $reads"
  goto done
endif

if ($mates == "") then
  set rr="-U $reads"
else
  set rr="-1 $reads -2 $mates"
endif
 
set type=`echo $reads | gawk '{t="f"}/fastq/{t="q"}{print t}'`

set params=""
if ($mm == 20_HISAT2_relaxed) then
  set params="--min-score L,0.0,-2"
endif

if (-e  $out.sam || -e $out.sam_sorted.gz) then
  echo "$out.sam ready"
else
  if (! -d $mm/$run) mkdir -p  $mm/$run
  uname -a
  echo " time bin/hisat2-master/hisat2 -$type -x $genomeDir/genome $rr $params -p $nThreads -S $out.sam"
	 time bin/hisat2-master/hisat2 -$type -x $genomeDir/genome $rr $params -p $nThreads -S $out.sam
endif
   
if (-e $out.sam && ! -e $out.sam_sorted.gz) then
  if (-e $TMPDIR/$mm/$run) \rm -rf $TMPDIR/$mm/$run
  mkdir -p $TMPDIR/$mm/$run
  cat   $out.sam | sort -T $TMPDIR/$mm/$run | gzip >  $out.sam_sorted.gz
endif
if (-e $out.sam && -e $out.sam_sorted.gz) then
   \rm  $out.sam
endif

if (-e $TMPDIR/$mm/$run) \rm -rf $TMPDIR/$mm/$run

done:
 echo done
