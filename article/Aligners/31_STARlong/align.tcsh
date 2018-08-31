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

set nThreads=4
set genomeDir=Aligners/$mm/$target
if (! -e $genomeDir) mkdir -p $genomeDir

if (! $?TMPDIR) then
  set TMPDIR=/tmp
  if (-d /export/home/TMP) set TMPDIR=/export/home/TMP
endif

# --sjdbOverhang 100 : best for long reads >= 100, do not give is no gff file
if (-e $genome  && ! -e $genomeDir/SAindex) then
  echo "gunzip -c $genome > $genomeDir/mygenome.fasta"
        gunzip -c $genome > $genomeDir/mygenome.fasta
  ls -ls $genomeDir/mygenome.fasta
  if (-e $genomeDir/mygenome.fasta) then
    # lmem03 8threads 1h49 elapsed, 513%, 32673u+1117s
    time bin/STARlong --runMode genomeGenerate --runThreadN $nThreads --genomeDir $genomeDir --genomeFastaFiles  $genomeDir/mygenome.fasta 
    touch $genomeDir/done
    #\rm  $genomeDir/mygenome.fasta
  endif
endif

if (-e $genomeDir/SAindex) then
  echo "STAR index is ready"
else
  echo "missing STAR index $genomeDir/SAindex"
  goto done
endif

if (0) then
--readFilesCommand gunzip -c # allow .gz on input
--outSAMtype SAM/BAM/None [Unsorted/SortedByCoordinates]
--outSAMattributes All
--outSAMattributes Standard
--outFileNamePrefix $out
--outTmpDir $TMPDIR
--outFilterMatchNmin 24
--outFilterScoreMin 24
--outFilterMismatchNmax 100000
--outFilterMismatchNoverLmax 100000
--twopassMode Basic
--genomeLoad LoadAndRemove
--genomeLoad NoSharedMemory  # only one compatible with 2-pass mode
endif
echo $out'Aligned.out.sam'

if (-e $out'Aligned.out.sam' && ! -e $out.sam_sorted.gz) then
  cat   $out'Aligned.out.sam' | sort -T $TMPDIR | gzip >  $mm/$run/$mm.$run.sam_sorted.gz
endif

set mySTAR=STAR
if ($mm == 31_STARlong) set mySTAR=STARlong
if ($mm == 30_STAR) then
  set mySTAR=STAR_1pass
endif

if (-e  $out.sam || -e $out.sam_sorted.gz) then
  echo "$out.sam ready"
else
  if (! -d $mm/$run) mkdir -p  $mm/$run
  mkdir -p $TMPDIR/$mm
  if (-e $TMPDIR/$mm/$run) \rm -rf $TMPDIR/$mm/$run
  if ($mm == 30_STAR) then
    uname -a
    echo " time bin/$mySTAR --runThreadN $nThreads --genomeDir $genomeDir --readFilesCommand gunzip -c --outSAMtype SAM --outSAMattributes Standard --outFileNamePrefix $out --outTmpDir $TMPDIR/$mm/$run --genomeLoad NoSharedMemory  --readFilesIn $reads $mates"
           time bin/$mySTAR --runThreadN $nThreads --genomeDir $genomeDir --readFilesCommand gunzip -c --outSAMtype SAM --outSAMattributes Standard --outFileNamePrefix $out --outTmpDir $TMPDIR/$mm/$run --genomeLoad NoSharedMemory  --readFilesIn $reads $mates
  else
    uname -a
    echo " time bin/$mySTAR --runThreadN $nThreads --genomeDir $genomeDir --readFilesCommand gunzip -c --outSAMtype SAM --outSAMattributes Standard --outFileNamePrefix $out --outTmpDir $TMPDIR/$mm/$run --outFilterMatchNmin 24 --outFilterScoreMin 24 --outFilterMismatchNmax 100000 --outFilterMismatchNoverLmax .5  --genomeLoad NoSharedMemory --twopassMode Basic --seedPerReadNmax 100000 --readFilesIn $reads $mates"
           time bin/$mySTAR --runThreadN $nThreads --genomeDir $genomeDir --readFilesCommand gunzip -c --outSAMtype SAM --outSAMattributes Standard --outFileNamePrefix $out --outTmpDir $TMPDIR/$mm/$run --outFilterMatchNmin 24 --outFilterScoreMin 24 --outFilterMismatchNmax 100000 --outFilterMismatchNoverLmax .5  --genomeLoad NoSharedMemory --twopassMode Basic --seedPerReadNmax 100000 --readFilesIn $reads $mates
  endif
endif
   

if (-e $out'Aligned.out.sam' && ! -e $mm/$run/$mm.$run.sam_sorted.gz) then
   cat   $out'Aligned.out.sam' | sort -T $TMPDIR | gzip >  $mm/$run/$mm.$run.sam_sorted.gz
endif
if (-e $out'Aligned.out.sam' && -e $mm/$run/$mm.$run.sam_sorted.gz) then
   \rm  $out'Aligned.out.sam' 
endif


done:
 echo done
