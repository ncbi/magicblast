#!/bin/tcsh

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
