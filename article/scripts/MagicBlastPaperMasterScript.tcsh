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
#
##       MagicBLAST_paper_master_script.tcsh
##
## MagicBLAST paper, june 2018, master script
## Author, Greg Boratyn, Danielle Thierry-Mieg, Jean Thierry-Mieg, Ben Busby, Tom Madden
## email for this script:   mieg@ncbi.nlm.nih.gov

## This is a tcsh executable script
## To see the on line help, run it under the tcsh interpretor using the command
##       MagicBLAST_paper_master_script.tcsh

if ($# == 0) goto phase_Help
if ($1 == help) goto phase_Help
if ($1 == '-help') goto phase_Help
if ($1 == '--help') goto phase_Help

#############################################################################
## Metadata

## Aligners
# List of aligners used in the analysis
# The number in front serves to order the tables in a systematic way
# one can insert new versions of each program by inserting new numbers
# but since the numbers are erased in the final tables, the number AND the names 
# must be unique
setenv methods "10_MagicBLAST  20_HISAT2_relaxed 21_HISAT2 30_STAR 31_STARlong  32_STAR.2.6c 40_TopHat2"
# setenv methods "21_HISAT2 31_STARlong"
# setenv methods "21_HISAT2"
# setenv methods "40_TopHat2"
# setenv methods "20_HISAT2_relaxed 21_HISAT2 "
# setenv methods "10_MagicBLAST"
# setenv methods "30_STAR 31_STARlong  32_STAR.2.6c"
# setenv methods "10_MagicBLAST  20_HISAT2_relaxed 21_HISAT2 30_STAR 31_STARlong  32_STAR.2.6c"
#############################################################################
## Datasets
## Each dataset must be aligned on the reference genome carrying the relevant truth
# Experimental human datasets, to be aligned on the NCBI human genome
setenv main_runs "iRefSeq PacBio Roche Illumina"

# Baruzzo datasets, to be aligned on the Baruzzo human reference genome
setenv HG19_r1_runs "HG19t1r1 HG19t2r1 HG19t3r1"
setenv HG19_r2_runs "HG19t1r2 HG19t2r2 HG19t3r2"
setenv HG19_r3_runs "HG19t1r3 HG19t2r3 HG19t3r3"
setenv HG19_runs "$HG19_r1_runs $HG19_r2_runs $HG19_r3_runs"

# Baruzzo datasets, to be aligned on the Baruzzo malaria reference genome
setenv PFAL_r1_runs "PFALt1r1 PFALt2r1 PFALt3r1"
setenv PFAL_r2_runs "PFALt1r2 PFALt2r2 PFALt3r2"
setenv PFAL_r3_runs "PFALt1r3 PFALt2r3 PFALt3r3"
setenv PFAL_runs "$PFAL_r1_runs $PFAL_r2_runs $PFAL_r3_runs"

setenv runs "$main_runs $HG19_runs  $PFAL_runs"

# Additional PacBio runs from Brain and Testis
setenv pacbio_runs "SRR5189652 SRR5189667"
# Additional long paired end Illumina reads 
# 250_250 (from metastatic melanoma)  and 300+300 from MCF7 cells)
setenv long_illumina_runs "SRR5438850 SRR5437876"
# setenv runs "$main_runs $HG19_runs  $PFAL_runs $pacbio_runs $long_illumina_runs"
# setenv runs "$long_illumina_runs"
# setenv runs "$pacbio_runs"
# setenv runs "PacBio Illumina"
# setenv runs "Roche PacBio iRefSeq "
# setenv runs "$runs PFALt1r1S HG19t1r1_50 $pacbio_runs"
# setenv runs "HG19t1r1_50"
# setenv runs "PFALt1r1S"
# setenv runs "$pacbio_runs"
# setenv runs "$runs HG19t1r1_50 PFALt1r1S"
setenv runs "$main_runs $HG19_runs  $PFAL_runs $pacbio_runs $long_illumina_runs PFALt1r1S HG19t1r1_50 "


# This adapter is present in the PacBio SRR runs and gives a peak at 32 aligned bases = polyA + first 8 bp of adaptor
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACTCT  GCGTTGATACCACTGCTTAGATCGGAAGAG
#############################################################################
## Fasta and Fastq files
## All runs fasta or fastq files are in the directories Fasta/$run
##  If they are absent, the script will download them from NCBI
#
#  The script assumes that all files are gzipped, and called $run/$run*.fast[aq].gz
#  their logical name, i.e. PacBio, links to their SRA identificator, e.g. SRR5009494.
#    The iRefSeq and the Baruzzo files are given in fasta format
#    The Illumina, Pabio and Roche fils are given in .fastq format
#    Some runs are paired-ends:
#       -Illumina paired end run has 2 files called SRR534301_1 and SRR534301_2
#       -Baruzzo paired end runs are called .forward and .reverse
#    In the Roche file Fasta/Roche/SRR899420.fastq we removed all read_1 (all 4-bases long)
#      and kept only the 24577 read_2.
#
#  For convenience and completeness, we also copied in Fasta/iRefSeq the original gff file
#  The iRefSeq fasta file can be extracted from the gff file using the command option
#       MagicBLAST_paper_master_script.tcsh make_iRefSeq 
#  of the the present script

foreach run ($runs)
  if (! -d Fasta/$run) mkdir -p  Fasta/$run
end 

#############################################################################
## Reference genome
# The main genome is NCBI release 19, limited to the main chromosome and 
# the mitochondria excluding the patches and the alternates
setenv main_genome GRCh38.genome.fasta.gz
# Baruzzo benchmark reference genome, available from
setenv HG19_genome HG19.Baruzzo.genome.fasta.gz
setenv PFAL_genome PFAL.Baruzzo.genome.fasta.gz

## Automatic download of the BenchMark fastq files from NCBI
set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/"
if (! -d Reference_genome) mkdir Reference_genome
pushd Reference_genome
  foreach ff (GRCh38.genome.fasta.gz  HG19.Baruzzo.genome.TM.txt  HG19.Baruzzo.genome.fasta.gz  PFAL.Baruzzo.genome.TM.txt  PFAL.Baruzzo.genome.fasta.gz)
    if (! -e $ff) then
      wget $FTP/Reference_genome/$ff
    endif
  end
popd

if (1) then
  foreach run ($HG19_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/HG19.Baruzzo.genome.fasta.gz genome.gz
      echo HG19 > target
    popd
  end
  foreach run ($PFAL_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/PFAL.Baruzzo.genome.fasta.gz genome.gz
      echo PFAL > target
    popd
  end
  foreach run ($main_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/GRCh38.genome.fasta.gz genome.gz
      echo GRCh38 > target
    popd
  end
  foreach run ($pacbio_runs $long_illumina_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/GRCh38.genome.fasta.gz genome.gz
      echo GRCh38 > target
    popd
  end
  touch Fasta/genomes_are_assigned

  foreach run ($pacbio_runs iRefSeq PacBio)
    if (-d  Fasta/$run) touch Fasta/$run/isLongRun
  end
endif

#############################################################################
## Automatic download of the binaries from the NCBI ftp site

set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article"
if (! -d bin || ! -d scripts/HTSeq) then
  if (-e bin/binaries.linux64.tar.gz) then
    mv binaries.linux64.tar.gz .
  endif
  if (! -e binaries.linux64.tar.gz) then  
     wget $FTP/binaries.linux64.tar.gz
  endif
  if (! -e binaries.linux64.tar.gz) then  
    echo "FATAL ERROR: The automatic download of the binaries from $FTP/binaries.linux64.tar.gz failed"
    echo "May be the network connection did not work, please try manually to run the command"
    echo "    wget $FTP/binaries.linux64.tar.gz"
    echo "if it does not work please email mieg@ncbi.nlm.nih.gov"
  endif
  if (-e binaries.linux64.tar.gz) then
    echo "expanding binaries.linux64.tar.gz, please wait"
    gunzip -c binaries.linux64.tar.gz | tar xf -
    mv binaries.linux64.tar.gz bin
  endif
endif

#############################################################################
## BAM files
## The BAM files are named    $method/$run/$method.$run.bam
## All runs were aligned on the relevant appropriate genome by all aligners
## but it did not always work. Some files are missing, for example 30_STAR.iRefSeq.bam, 
## because STAR crashed on long reads

##############################################################################
## utilities
setenv TMPDIR /tmp
if (-d /export/home/TMP) setenv TMPDIR /export/home/TMP
if (! -d tmp) mkdir tmp
if (! -d RESULTS) mkdir RESULTS

##############################################################################
##############################################################################
## Executable and source code
## Our scripts are in the scripts directory, e.g. scripts/AliQC.py
## Our executable are compiled for generic LINUX 64 bits machine in the bin directory
##    e.g. magicblast, dna2dna, sam2gold 
## Our source code is available for analysis and recompilation in machine optimized mode
## in the source_code directory, together with instructions in the correspondng README file.

echo -n "## MagicBlastPaperMasterScript.tcsh $1 : "
date

echo "runs = $runs"
echo "methods = $methods"

if ($1 == init) goto phase_Init
if ($1 == download) goto phase_Download
if ($1 == align) goto phase_Align
if ($1 == Make_iRefSeq) goto phase_Make_iRefSeq
if ($1 == count) goto phase_Count
if ($1 == sam) goto phase_Sam
if ($1 == accuracy) goto phase_Accuracy
if ($1 == aliqc) goto phase_aliqc
if ($1 == errors) goto phase_DirectErrorCount
if ($1 == export) goto phase_Export
if ($1 == aliLn) goto phase_aliLn
if ($1 == subs) goto phase_Count_subtitutions_in_benchmark
echo "Unknown command : $1, please try --help"
goto phaseLoop

phase_Help:

echo "\nusage  scripts/MagicBlastPaperMasterScript.tcsh, where command is one of:"
echo 'help   : This cnline help'
echo 'init   : in tcsh, "source README init" will set the variables $runs, $methods which may be convenient'
echo 'download  : Automatic download of the fastq files, please monitor carefully the results'
echo 'Make_iRefSeq : create the iRefSeq fasta file and intron file from the gff and the genome file'
echo 'count  : Count the reads in each fasta/fastq file'
echo 'sam    : Automatically download the sam files from NCBI (rather that running the aligners locally)'
echo 'align  : run for all runs all aligners for which the script Aligners/$method/align.tcsh is defined'
echo 'aliqc : run AliQC.py on each run in the background (slow), presupposes that HTSeq is installed'
echo 'accuracy : measure intron and alignment accuracy relative to the gold standard'
echo 'errors : count the errors in te BAM files using the NM:i:x optional field'
echo 'subs   : count substitutions in the human and malaria Baruzzo benchmarks'
echo 'aliLn  : export histogram of aligned lengths'
echo 'export : export QC and ROC curve of intron discovery'

goto phaseLoop

phase_Init:
goto phaseLoop

##############################################################################
#############################################################################
## Automatic download of the BenchMark fastq files from NCBI
phase_Download:

## NCBI repository for the datafiles used in the paper
set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/"

echo "checking $HG19_runs  $PFAL_runs"
foreach run ( $HG19_runs  $PFAL_runs)
  if (! -e Fasta/$run/$run.reverse.fasta.gz) then
    pushd Fasta/$run
      echo "Trying to download $run from$FTP"
      wget $FTP/Fasta/$run/$run.cig.gz
      wget $FTP/Fasta/$run/$run.forward.fasta.gz
      wget $FTP/Fasta/$run/$run.reverse.fasta.gz
    popd
  endif
end

set run=HG19t1r1_50
if (! -d Fasta/$run) then
  echo "preparing the 5+50 clipped run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    ln -s ../../Reference_genome/HG19.Baruzzo.genome.fasta.gz genome.gz
    ../../bin/dna2dna -i ../HG19t1r1/HG19t1r1.forward.fasta.gz -gzo -o $run.forward -rightClipAt 50
    ../../bin/dna2dna -i ../HG19t1r1/HG19t1r1.reverse.fasta.gz -gzo -o $run.reverse -rightClipAt 50
  popd
endif

set run=PFALt1r1S
if (! -d Fasta/$run) then
  echo "preparing the subsampled run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    ln -s ../../Reference_genome/PFAL.Baruzzo.genome.fasta.gz genome.gz
    ../../bin/dna2dna -i ../PFALt1r1/PFALt1r1.forward.fasta.gz -gzo -o $run.forward -subsample 100
    ../../bin/dna2dna -i ../PFALt1r1/PFALt1r1.reverse.fasta.gz -gzo -o $run.reverse -subsample 100
  popd
endif

echo "checking iRefSeq"
foreach run (iRefSeq)
  if (! -e Fasta/$run/$run.fasta.gz) then
    pushd Fasta/$run
      echo "Trying to download $run from$FTP"
      wget $FTP/Fasta/$run/$run.cig.gz
      wget $FTP/Fasta/$run/$run.fasta.gz
      wget $FTP/Fasta/$run/GRCh38_genome.gff.gz
      ln -s GRCh38_genome.gff.gz genome.gff.gz
      gunzip -c $run.fasta.gz | ../../bin/dna2dna -getTM > $run.TM
    popd
  endif
end

echo "checking Roche"
foreach run (Roche)
  if (! -e Fasta/$run/$run.fasta.gz) then
    pushd Fasta/$run
      echo "Trying to download $run from$FTP"
      wget $FTP/Fasta/$run/$run.fasta.gz
      gunzip -c $run.fasta.gz | ../../bin/dna2dna -getTM > $run.TM
    popd
  endif
end

#############################################################################
## .cig TRUTH Files
## The Baruzzo benchmark is providing the original position of the simulated reads
## in their .cig format, which is analogous, but not identical, to the SAM format.
## Since the fasta and the .cig files both come from Baruzzo, we located then in Fasta/$run.cig.gz
## For convenience, we reformatted the RefSeq gff file into a similar Fasta/iRefSeq/iRefSeq.cig.gz
if (-e Fasta/iRefSeq/genome.gff.gz && ! -e Fasta/iRefSeq/iRefSeq.cig.gz) then
  gunzip -c Fasta/iRefSeq/genome.gff.gz | gawk -F '\t' -f scripts/gff2cig.awk | gzip >  Fasta/iRefSeq/iRefSeq.cig.gz4
endif
## To compare the BAM files produced by the different aligners to the .cig "truth"
## we developped a C code called sam2gold (see below)

#############################################################################
## Automatic download of the Illumina Roche pacBio from SRA

foreach run (PacBio Illumina)
  if (! -d Fasta/$run) continue
  if ($run == PacBio) set run2=SRR5009494
  if ($run == Roche)  set run2=SRR899420
  if ($run == Illumina)  set run2=SRR534301
  if (! -e Fasta/$run/$run2.fastq.gz && ! -e Fasta/$run/$run2'_1'.fastq.gz) then
     set n=`bin/fastq-dump --help | wc -l`
    if ($n < 10) then
      echo "Sorry, the executable bin/fastq-dump available from NCBI SRA and needed to download the $run run is not found"
      goto phaseLoop
    endif
    set sf=""
    if ($run == Illumina)  set sf="--split-files"
    echo "Trying to download $run2 from SRA"
    bin/fastq-dump $sf -O Fasta/$run $run2
    if (-e Fasta/$run/$run2.fastq || -e Fasta/$run/$run2'_1'.fastq) then
      gzip Fasta/$run/$run2*.fastq
      pushd Fasta/$run
        if (-e $run2.fastq.gz) ln -s $run2.fastq.gz $run.fastq.gz
        if (-e $run2'_1'.fastq.gz) ln -s $run2'_1'.fastq.gz $run'_1'.fastq.gz
        if (-e $run2'_2'.fastq.gz) ln -s $run2'_2'.fastq.gz $run'_2'.fastq.gz  
      popd
    endif
    if (-e ~/ncbi/public/sra/$run2) \rm -rf  ~/ncbi/public/sra/$run2
  endif
end

#############################################################################
## Automatic download of the fastq files from SRA

foreach run ($pacbio_runs)
  if (! -d Fasta/$run) continue
  if (! -e Fasta/$run/$run.fasta.gz && ! -e Fasta/$run/$run.fastq.gz) then
    set n=`bin/fastq-dump --help | wc -l`
    if ($n < 10) then
      echo "Sorry, the executable bin/fastq-dump available from NCBI SRA and needed to download the $run run is not found"
      goto phaseLoop
    endif
    echo "Trying to download $run from SRA"
    bin/fastq-dump -O Fasta/$run $run
    gzip Fasta/$run/$run.fastq
  endif
end

foreach run ($long_illumina_runs)
  if (! -d Fasta/$run) continue
  if (! -e Fasta/$run/$run.fasta'_1'.gz && ! -e Fasta/$run/$run'_1'.fastq.gz) then
    set n=`bin/fastq-dump --help | wc -l`
    if ($n < 10) then
      echo "Sorry, the executable bin/fastq-dump available from NCBI SRA and needed to download the $run run is not found"
      goto phaseLoop
    endif
    echo "Trying to download $run2 from SRA"
    bin/fastq-dump -O Fasta/$run --split-files $run 
    gzip Fasta/$run/$run*.fastq
  endif
end

goto phaseLoop

##############################################################################
##############################################################################
## Count the number of reads, the shortest, the longest read in every fasta/fastq file
## using the utility bin/dna2dna (compiled for Linux 64bits) 
## The source code is part of the aceview/magic distribution in the source_code directory

phase_Count:
echo 'counting the number of reads in each fasta/fastq file'
foreach run ($runs)
  if (-e Fasta/$run/$run.fasta.gz && ! -e Fasta/$run/$run.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run.fasta.gz -I fasta -count -o Fasta/$run/$run
  endif
  if (-e Fasta/$run/$run.fastq.gz && ! -e Fasta/$run/$run.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run.fastq.gz -I fastq -count -o Fasta/$run/$run
  endif
  if (-e Fasta/$run/$run'_1'.fastq.gz && ! -e Fasta/$run/$run'_1'.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run'_1'.fastq.gz -I fastq -count -o Fasta/$run/$run'_1'
    bin/dna2dna -i  Fasta/$run/$run'_2'.fastq.gz -I fastq -count -o Fasta/$run/$run'_2'
  endif
  if (-e Fasta/$run/$run.forward.fasta.gz && ! -e Fasta/$run/$run.forward.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run.forward.fasta.gz -I fasta -count -o Fasta/$run/$run.forward 
    bin/dna2dna -i  Fasta/$run/$run.reverse.fasta.gz -I fasta -count -o Fasta/$run/$run.reverse 
  endif
  set nreads=`cat Fasta/$run/*.count  | gawk '/^Fragment_kept/{n+=$2}END{print n}'`
  echo "$run contains $nreads reads"
end
 
goto phaseLoop

##############################################################################
##############################################################################
## Create the iRefSeq fasta file and intron file from the gff and the genome file

phase_Make_iRefSeq:

## In practice, the file Fasta/iRefSeq/iRefSeq.fasta.gz is downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/
## The scrip is given here for transparency and to allow the reconstruction
## of the iRefSeq in the future from a diferent gff file and reference genome

if (! -e Fasta/iRefSeq/iRefSeq.fasta.gz) then
  echo "Creating Fasta/iRefSeq/iRefSeq.fasta.gz using the genome and the gff3 anntation"
  if (! -e Fasta/iRefSeq/genome.gz) then
    echo "Missing file Fasta/iRefSeq/genome.gz, I cannot create the iRefSeq fasta file"
    goto phaseLoop
  endif
  if (! -e Fasta/iRefSeq/genome.gff.gz( then
     echo "Missing file Fasta/iRefSeq/genome.gff.gz, I cannot create the iRefSeq fasta file"
    goto phaseLoop
  endif

  echo "Found the genome and the gff file, constructing the fasta in Fasta/iRefSeq/tmp"
  if (! -d Fasta/iRefSeq/tmp) mkdir  Fasta/iRefSeq/tmp
  pushd  Fasta/iRefSeq/tmp
    # This script is surprisingly complex, sorry, because we are trying to identify the NMs which map as well
    # at different locus of the genome, but while doing so, we unfortunately discovered a number of
    # irregularities in the definition of the RefSeqs that we try to compensate

    # To simplify the matter, we directly provide the iRefSeq fasta and gff files on our ftp site. 

    # extract the NM_ from the gff file
    zcat ../genome.gff.gz | grep NM_ | grep NC_ > NM.gff
    # we could directly export the fasta file with the command 
    # ../../../bin/dna2dna -gff3 NM.gff -gtfRemap iRefSeq -gtfGenome ../genome.gz -o iRefSeq -O fasta
    # but some NM have a single indentifier and yet map on 2 chroms
    # by not providing the genome we only export the 6 columns sponge file
    ../../../bin/dna2dna -gff3 NM.gff -gtfRemap iRefSeq -o iRefSeq -O fasta
    set nNM=`cat iRefSeq.[fr].sponge | cut -f 1 | sort -u | wc -l`
    echo  "Number of NM_ $nNM (is 45065)"
    set nNM_chr=`cat iRefSeq.[fr].sponge | cut -f 1,3 | sort -u | wc -l`
    echo  "Number of NM_chrom $nNM_chr (is 45108)"
    set nG=`cat iRefSeq.[fr].sponge | cut -f 6 | sort -u | wc -l`
    echo  "Number of genes with NM $nG (is  19269)"
    echo "Evaluating the mapping multiplicity of the iRefSeq"
    # to fix the issue that the same NM may map on  several chromosomes 45108 = 45065 = 43 cases
    # we merge the chrom and the NM in column 1 to create a disambiguated shadow file
    cat  iRefSeq.[fr].sponge  | gawk -F '\t' '{printf("%s:%s\t%s\t%s\t%s\t%s\t%s\n",$1,$3,$2,$3,$4,$5,$6);}' > NM_chr.sponge
    # the sponge file has the the NM the gene and the coordinates of all exons, hence the sequence
    # we now construct the fasta file
     ../../../bin/dna2dna -shadow NM_chr.sponge -i ../genome.gz -o iRefSeq -O fasta -maxLineLn 80 -gzo
    # measure the number of disntinct NM with identical sponge (hence DNA) and mapping    
    cat iRefSeq.[fr].sponge | grep NM_ | sort > _t
    cat _t | gawk '{nm=$1;z[nm]=z[nm] "_" $3 "_" $4 "_" $5;}END{for(k in z){g[z[k]]=k;u[z[k]]+=1;}for (v in u) {if(u[v]>1)print u[v],g[v]}}' | sort -k 1n > iRefSeq.mapping_multiplicity
    ../../../dna2dna -i iRefSeq.fasta.gz  -count -o iRefSeq
    \mv iRefSeq.fasta.gz iRefSeq.count ..
    # extract
    ../../../dna2dna -i iRefSeq.fasta.gz -O raw | sort > _t1
    # count the distnct NM sequences : 44914
    ../../../dna2dna -i iRefSeq.fasta.gz -O raw | cut -f 1 | sort -u | wc
    ../../../dna2dna  -i iRefSeq.fasta.gz -getTM > ../iRefSeq.TM

    # map the NM on the NM to find who is identical or included in the other
    clipalign -i iRefSeq.fasta.gz -t iRefSeq.fasta.gz -errMax 0 -o nm2nm -maxHit 24 -minAli 140
    bestali -i  nm2nm.hits -exportBest -o nm2nm33
    # now count NM mapping exactly in NM with a different geneid -> 143, we add 43+43 for the 43 NM which map on X and Y
    cat nm2nm33.hits | gawk '{if($2-$4==0 && index($1,"|"$9">")==0)print}' > nm2nm.2genes.hits
    wc nm2nm.2genes.hits
    cat nm2nm.2genes.hits | gawk  -F '\t' '{n[$1]++;}END{for(k in n)u[n[k]]++;for (v  in u) {if(v>0)k+=u[v];kk+=v*u[v];print v, u[v];}print k, kk}' | sort -k 1n
    # we now have 303 NM mapping on another NM with a different gene name, however some distinct genes have same gene coordinates
    # extract the extreme coords of the NM from the sponge file
    cat NM_chr.sponge | gawk -F '\t' '{nm=$1;a1=$4;if($5<a1)a1=$5;a2=$5;if($4>a2)a2=$4;if(aa2[nm]<a2)aa2[nm]=a2;if(aa1[nm]>-a1)aa1[nm]=-a1;}END{for(nm in aa1) printf("%s\t%d\t%d\n",nm,-aa1[nm],aa2[nm]);}' > NM_chr.segment
    # reanalize the nm 2 nm hits file and eliminate the lines with overlapping coordinates
    echo ZZZZZ > ZZZZZ
    cat NM_chr.segment ZZZZZ  nm2nm.2genes.hits | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{nm=$1;if(zz+0<1){aa1[nm]=$2;aa2[nm]=$3;split(nm,aa,":");chrom[nm]=aa[2];next;}}{split($1,aa,"|");nm1=aa[1];nm2=$11;ok=1;if (chrom[nm1]==chrom[nm2] && aa1[nm1]<aa2[nm2] && aa2[nm1] > aa1[nm2])ok=0;if (ok==1)print}' > nm2nm.2genes.hits.no_doublons
    cat NM_chr.segment ZZZZZ  nm2nm.2genes.hits | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{nm=$1;if(zz+0<1){aa1[nm]=$2;aa2[nm]=$3;split(nm,aa,":");chrom[nm]=aa[2];next;}}{split($1,aa,"|");nm1=aa[1];nm2=$11;ok=1;if (chrom[nm1]==chrom[nm2] && aa1[nm1]<aa2[nm2] && aa2[nm1] > aa1[nm2])ok=0;if (ok==0)print}' > nm2nm.2genes.hits.doublons
    # final count of the repeated NM : 291 NM have several mappings
    cat nm2nm.2genes.hits.no_doublons | gawk  -F '\t' '{n[$1]++;}END{for(k in n)u[n[k]]++;for (v  in u) {if(v>0)k+=u[v];kk+=v*u[v];print v, u[v];}print kk}' | sort -k 1n
    # so finally we have 291 NM have multiple mapping just by looking at the annotated NM themselves + (43 + 43)  from the pseudo autosomal region with single NM and geneid total 291+86=379
    cat NM_chr.segment| gawk '{split($1,aa,":");n[aa[1]]++;chrom[aa[1]]=aa[2];}END{for (nm in n)if(n[nm]>1)print nm,n[nm],chrom[nm];}' > NM.pseudo_autosomal_region.mapping_twice 
    wc  NM.pseudo_autosomal_region.mapping_twice     
    cat nm2nm.2genes.hits.no_doublons | gawk '{split($1,aa,"|");print aa[3] "="$9}' | sed -e 's/>/</' | sort -u > gene_pairs
    ## construct a cig file for the refseq
    # Use NM_...:chrom... as NM identifiers because in the pseudo autosomal region, the same NM maps in 2 places: one NM_ 2 locus
    # whereas a usual palindromic exactly duplicated genes has 1 NM per locus, i.e. 2 NM 2 locus
    # this raises the number of NM supported introns from 210357 to 210509
    
    zcat ../genome.gff.gz | gawk -F '\t' '{if ($3 != "exon") next;}{split($9,aa,"Genbank:");split(aa[2],bb,",");split(bb[1],cc,";");seq=cc[1];if(substr(seq,1,2)!="NM" && substr(seq,1,2)!="zNR")next;seq=seq ":" $1; chrom[seq]=$1;nx[seq]++;i=nx[seq];a1[seq,i]=$4;a2[seq,i]=$5;strand[seq]=$7;}END{for (seq in nx){n=nx[seq];printf("%s\t%s",seq,chrom[seq]);if(strand[seq]=="+"){printf("\t%d\t%d\t", a1[seq,1], a2[seq,n]) ;for(i = 1 ; i <=n ; i++){if(i>1){dx=a1[seq,i]-a2[seq,i-1]-1;printf("%dN",dx);}dx=a2[seq,i]-a1[seq,i]+1;printf("%dM",dx);}}else {printf("\t%d\t%d\t", a1[seq,n], a2[seq,1]) ;for(i = n ; i >=1 ; i--){if(i<n){dx=a1[seq,i]-a2[seq,i+1]-1;printf("%dN",dx);}dx=a2[seq,i]-a1[seq,i]+1;printf("%dM",dx);}}printf("\t.\t%s\t.\n",strand[seq]);}}' | gzip > iRefSeq.cig.gz
    mv iRefSeq.cig.gz ..

endif

goto phaseLoop

##############################################################################
##############################################################################
## SAM 
## download the precomputed SAM files from NCBI
phase_Sam:

set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/"
foreach run ($runs)
  foreach method ($methods)

    ## the preferred methos is to download the aligned files from NCBI

    # For HISAT and STAR we have a special version of the code to align long runs
    # so in these cases we do not atempt to align the long runs with the short code
    if ($method == 30_STAR || $method == 32_STAR.2.6c) then
       if (-e Fasta/$run/isLongRun) continue
    endif
    # and vice versa
    if ( $method == 20_HISAT2_relaxed) then
       if (! -e Fasta/$run/isLongRun) continue
    endif
 
      if (! -e $method/$run/$method.$run.sam_sorted.gz) then
        mkdir -p $method/$run
        pushd  $method/$run
          wget $FTP/SAM/$method.$run.sam_sorted.gz
        popd
      endif

  end
end

goto phaseLoop

##############################################################################
##############################################################################
## ALIGN Run all aligners on all runs

phase_Align:

foreach run ($runs)
  foreach method ($methods)

    if (! -e Aligners/$method/align.tcsh) then
      echo "missing script Aligners/$method/align.tcsh"
      continue
    endif

    # For HISAT and STAR we have a special version of the code to align long runs
    # so in these cases we do not atempt to align the long runs with the short code
    if ($method == 30_STAR || $method == 32_STAR.2.6c) then
       if (-e Fasta/$run/isLongRun) continue
    endif
    # and vice versa
    if ( $method == 20_HISAT2_relaxed) then
       if (! -e Fasta/$run/isLongRun) continue
    endif

    if (-e Aligners/$method/align.tcsh && ! -e $method/$run/$method.$run.sam_sorted.gz) then
      set read_1="x"
      set read_2=""
      if (-e Fasta/$run/$run'_1.fastq.gz') set read_1=Fasta/$run/$run'_1.fastq.gz'
      if (-e Fasta/$run/$run'_2.fastq.gz') set read_2=Fasta/$run/$run'_2.fastq.gz'
      if (-e Fasta/$run/$run'_1.fasta.gz') set read_1=Fasta/$run/$run'_1.fasta.gz'
      if (-e Fasta/$run/$run'_2.fasta.gz') set read_2=Fasta/$run/$run'_2.fasta.gz'
      if (-e Fasta/$run/$run.fasta.gz) set read_1=Fasta/$run/$run.fasta.gz
      if (-e Fasta/$run/$run.fastq.gz) set read_1=Fasta/$run/$run.fastq.gz
      if (-e Fasta/$run/$run.forward.fastq.gz) set read_1=Fasta/$run/$run.forward.fastq.gz
      if (-e Fasta/$run/$run.reverse.fastq.gz) set read_2=Fasta/$run/$run.reverse.fastq.gz
      if (-e Fasta/$run/$run.forward.fasta.gz) set read_1=Fasta/$run/$run.forward.fasta.gz
      if (-e Fasta/$run/$run.reverse.fasta.gz) set read_2=Fasta/$run/$run.reverse.fasta.gz

      if (! -e $read_1) then
        echo "Run $run Missing read file $read_1"
	ls -ls Fasta/$run/*fast*
        continue
      endif
      set target=`cat Fasta/$run/target`
      if (! -d $method/$run) mkdir -p $method/$run
        echo "align $method/$run"
        scripts/submit $method/$run "Aligners/$method/align.tcsh $method $run $target $read_1 $read_2" 
        # scripts/submit $method/$run "Aligners/$method/align.tcsh $method $run $target $read_1 $read_2" 64G UGE4
    endif
  end
end

goto phaseLoop

##############################################################################
##############################################################################
## Intron, exon, insertion deletion comparison to the TRUTH
## Compare the alignment results, provided in BAM format 
## to the GOLD standard truth from iRefSeq and Baruzzo given in .cig format
## The source C-code is part of the aceview/magic distribution www.aceview.org/Software
## The executable for LINUX 64 bits is in bin
##
## sam2gold produces several output files
##    .qc a small self documented statistics table in tab delimited format
##    .aliqc the same statistics in a more computer friendly tag-values tab delimited format
##    .Intron a table giving the coordinates of all introns, with support in GOLD or BAM
##    .Deletion a table giving the coordinates of all deletions, with support in GOLD or BAM
##    .Insertion a table giving the coordinates of all insertions, with support in GOLD or BAM

phase_Accuracy:
foreach run ($runs)
  if (-e Fasta/iRefSeq/iRefSeq.cig.gz && ! -e Fasta/$run/$run.cig.gz) then
    pushd Fasta/$run
      ln -s ../iRefSeq/iRefSeq.cig.gz   $run.cig.gz
    popd
  endif
end
# Illumina $HG19_runs  $PFAL_runs 
#$main_runs $pacbio_runs $long_illumina_runs  $HG19_runs  $PFAL_runs $methods
foreach run ($runs)
  foreach mm ($methods)
     if (-e Fasta/$run/$run.cig.gz && -e $mm/$run/$mm.$run.sam_sorted.gz && ! -e $mm/$run/$mm.$run.delins.tsv) then
       echo $mm/$run/$mm.$run.sam_sorted.gz
       set arp=""
       set arp=`echo $mm | gawk 'BEGIN{arp="";}{if(index($1,"STAR")>0) arp="-addReadPairSuffix"}END{printf("%s",arp);}'`
       \rm  $mm/$run/$mm.$run.sam2gold.*
       scripts/submit $mm/$run/$mm.$run.sam2gold "bin/sam2gold $arp -g $run..GOLD:Fasta/$run/$run.cig.gz -i  $run..$mm':'$mm/$run/$mm.$run.sam_sorted.gz -o  $mm/$run/$mm.$run"
      endif
  end
end
 
goto phaseLoop

##############################################################################
##############################################################################
## Alignment quality control
## Evaluate in great details the intrinsic quality of the alignment results, provided in BAM format 
## This analysis does not refer to the gold truth
## This is a python.2.7 scripts given in scripts/AliQC.py
## It was developped in collaboration with Joe Meehan, FDA, for the MAQC/SEQC project
## There is a dependency, one must fisrt install HTSeq as explained in the previous section
##
## aliqc produces a computer friendly tag-values tab delimited format .aliqc.tsv
## In the following section aliqc is used again to merge these file into a single table

phase_aliqc:

# create sam_sorted only once
set ok=1
  foreach mm ($methods)
    foreach run ($runs)
      if (-e $mm/$run/$mm.$run.bam && ! -e $mm/$run/$mm.$run.sam_sorted.gz) then
        set ok=0
        echo "transformng $mm/$run/$mm.$run.bam into .sam_sorted.gz"
        scripts/submit $mm/$run/$mm.$run.samview "samtools view $mm/$run/$mm.$run.bam | sort -T $TMPDIR | gzip >  $mm/$run/$mm.$run.sam_sorted.gz" 
      endif
    end
  end
if ($ok == 0) goto phaseLoop

# use sam_sorted.gz rather than bam (aliqc --BAM also works, but it need to call sort which is very costly)
# PacBio Roche iRefSeq SRR5189652 SRR5189667 HG19t1r1 HG19t1r2 HG19t2r1 HG19t2r2 HG19t3r1 PFALt1r1  PFALt1r2  PFALt1r3  PFALt2r1 PFALt2r2 PFALt2r3 PFALt3r1 PFALt3r2 PFALt3r3 SRR5438850
foreach minAli (0 50 80)
  foreach run ($runs)
    foreach mm ($methods)
      if (-e $mm/$run/$mm.$run.sam_sorted.gz && -e Fasta/$run/genome.gz && ! -e $mm/$run/$mm.$run.minAli$minAli.aliqc.tsv) then
        echo "Running AliQC.py on $mm/$run/$mm.$run.sam_sorted.gz"
	set nreads=`cat Fasta/$run/*.count  | gawk '/^Fragment_kept/{n+=$2}END{print n}'`
	if ($nreads == 0) then
	  echo "missing file Fasta/$run/$run.count, please run phase 1"
        else
          echo "Running AliQC.py $mm $run"
          scripts/submit $mm/$run/aliqc.minali$minAli "python3 scripts/AliQC.py --SAMSORTEDGZ -i $mm/$run/$mm.$run.sam_sorted.gz -r $run..$mm.minAli$minAli -f Fasta/$run/genome.gz -o $mm/$run/$mm.$run.minAli$minAli --nreads $nreads --minAli $minAli"
         endif
      endif
    end
  end
end

goto phaseLoop

##############################################################################
##############################################################################
## Direct statistics of the error counts reported in the BAM files
## The AliQC.py code, above, parses the bam file and the genome
## and computed its own evaluation of the number of error per aligned read
## In the present script, we rely on the presence in the BAM files of the NM:i:x
## optional field, collate the X values and report the statistics
## Hopefully, the 2 methods should be compatible, but they do not have
## to agree exactly since aliqc counts a double or triple deletion as a single event
## and some aligners may report it as 2 or 3 errors

phase_DirectErrorCount:
echo phase_DirectErrorCount
foreach mm ($methods)
  foreach run ($runs)
    echo "phase_DirectErrorCount $mm $run"
    if (-e $mm/$run/$mm.$run.sam_sorted.gz && ! -e $mm/$run/$mm.$run.nerrors) then
      scripts/submit $mm/$run/$mm.$run.nerrors "scripts/directErrorCount.tcsh $mm $run" 
     endif
  end
end

# export the results in a single human readable table
set toto=RESULTS/Error_counts_using_NM_optional_tag.txt
echo -n "## $toto :" > $toto
if (-e $toto.1) \rm $toto.1
date >> $toto
foreach mm ($methods)
  foreach run ($runs)
    if (-e  $mm/$run/$mm.$run.nerrors) then
      echo -n "$mm\t$run\t" >> $toto.1
      cat $mm/$run/$mm.$run.nerrors | gawk '/^#/{next;}{n[$1]=$2;if($1>max)max=$1;}END{for(i=-5;i<=max;i++)printf("\t%d",n[i]);printf("\n")}' >> $toto.1
    endif
  end
end
cat $toto.1 | gawk '{if(NF>max)max=NF;}END{printf("\t\t");max-=6;for(i=-5;i<=max;i++)printf("\t%d",i);printf("\n")}' > $toto.2
cat $toto.2 $toto.1 | gawk -F '\t' -f scripts/transpose.awk >> $toto
\rm $toto.1 $toto.2
echo "The table of errors using the optional NM tag of the BAM files is in"
echo $toto

goto phaseLoop

##############################################################################
##############################################################################
## Count the substitutions as declared in the Baruzzo Benchmark
## The statistics only measures the substitutions in the r3 reads
## fully and continuously aligned on the plus strand of chromosome 8
## This seems sufficient since it involves in each case at least 100,000 reads 

phase_Count_subtitutions_in_benchmark:

# there are several phases in the calsulation
# 1: select the full reads, characterized by a cigar string 100M
if (! -e SUBS) mkdir SUBS
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    if (-e Fasta/$sp$tt'r3'/$sp$tt'r3'.cig.gz && ! -e  SUBS/subs.$sp$tt) then 
      zcat Fasta/$sp$tt'r3'/$sp$tt'r3'.cig.gz | grep chr8  | grep '+' |  grep 100M | cut -f 1,2,3,4,8 > SUBS/subs.$sp$tt
    endif
  end
end

# 2: construct a 6 columns tab delimited shadow file, summarizing the coordinate of the alignemnts
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    if (-e  SUBS/subs.$sp$tt && ! -e SUBS/subs.$sp$tt.shadow ) then 
      cat SUBS/subs.$sp$tt | gawk -F '\t' '{printf("%s\t1\t100\t%s\t%d\t%d\n",$1,$2,$3,$4);}' > SUBS/subs.$sp$tt.shadow 
    endif
  end
end

# 3: isolate the genome of chromosome 8, using the dna2dna utility
foreach sp (HG19 PFAL)
  if (-e Reference_genome/$sp.Baruzzo.genome.fasta.gz && ! -e Reference_genome/$sp.Baruzzo.chr8.fasta.gz) then
    bin/dna2dna -i Reference_genome/$sp.Baruzzo.genome.fasta.gz -I fasta -O fasta -keepName -o Reference_genome/$sp.Baruzzo.chr8 -gzo 
  endif
end

# 4: Export the corresponding genomic segement in raw format
# The raw format has just 2 tab delimited columns: atgcatgc  identifier
# Notice that dna2dna is very versatile, it can directly export messenger RNAs given a genome and a gff file.
# try bin/dna2dna --help for a full list of functioalities
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    if (-e  SUBS/subs.$sp$tt.shadow && -e Reference_genome/$sp.Baruzzo.chr8.fasta.gz && ! -e SUBS/subs.$sp$tt.raw) then
      dna2dna -i Reference_genome/$sp.Baruzzo.chr8.fasta.gz -shadow SUBS/subs.$sp$tt.shadow -O raw -keepName >  SUBS/subs.$sp$tt.raw 
    endif
  end
end

# 5: the first subs file contains in column 1 and 5 the name and sequence of each read
#    the raw file contains in column 2 and 1 the name and sequence of each corresponding genomic segment
#    Both sequences are exactly 100 bases long, hence a simple awk script is sufficient to count the mismatching bases
echo ZZZZZ > SUBS/ZZZZZ
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    cat SUBS/subs.$sp$tt SUBS/ZZZZZ  SUBS/subs.$sp$tt.raw | gawk -F '\t' '/^ZZZZZ/{zz=1;}{if(zz<1){seq[$1]=$5;next;}if (seq[$2]){s1=seq[$2];s2=toupper($1);n=0;for(i=1;i<=100;i++)if(substr(s1,i,1) != substr(s2,i,1))n++;print n}}' | tags | sort -k 1n > SUBS/subs.$sp$tt.txt &
  end
end

# 6: produce a final synthetic table
set toto=RESULTS/mm_stats.Baruzzo.txt
echo -n "### $toto : " > $toto
date >> $toto
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    cat  SUBS/subs.$sp$tt.txt | gawk '{n1+=$2;n2+=$1*$2;printf ("%s\t%d\t%d\t%d\t%d\n",t,$1,$2,n1,n2);}' t=$sp$tt >> $toto
  end
end
echo "The distribution of substitutions in chromosome 8 in tbe Baruzzo benchmark"
echo "have been exported in $toto"
head -12 $toto

goto phaseLoop

##############################################################################
##############################################################################
## Exportation of global, human readable, quality control tables
## These tables where used directly to prepare the plots and table of the
## Magic-BLAST paper

phase_Export:

if (! -d RESULTS) mkdir RESULTS

# collate the aliqc.tsv tables from all runs using again the AliQC.py scripts with --table option
# this will produce 3 output tables
foreach minAli (0 50 80)
  if (-e toto) \rm toto
  foreach method ($methods)
    cat $method/*/*.minAli$minAli.aliqc.tsv  >> toto
  end
  cat toto | python3 scripts/AliQC.py --view table --split --minAli $minAli -o RESULTS/magicblast.table.minAli$minAli
end

# reformat the 3 output tables
foreach minAli (0 50 80)
  foreach type (mismatch_histo_and_types mismatch_per_cycle_per_kb_aligned aligned_reads_multiplicity_length_aligned)
    set toto=RESULTS/magicblast.table.minAli$minAli.$type
    cat $toto.tsv  | head -12 |  gawk -F '\t' '/^##/{next;}/^#/{print}' > $toto.txt
    cat $toto.tsv  | gawk -F '\t' '/^###/{next;}/^##/{print}' >> $toto.txt
    cat $toto.tsv  | gawk -F '\t' '/^#/{next;}{print}' | sort -k 1,1 -k 2,2 -k 3,3 | grep -v r2 | grep -v r3 | sed -e 's/r1//g' |  gawk -F '\t' '{if ($2 != old2){printf("\n");old2=$2;old3="";}}{if ($3 != old3){printf("\n");old3=$3}}{print}' >> $toto.txt
    cat $toto.tsv  | gawk -F '\t' '/^###/{print}' >> $toto.txt
  end
end

set toto=RESULTS/Mapping_accuracy.txt
echo -n "## $toto :" > $toto
date >> $toto 
cat scripts/mapping_accuracy.header   >> $toto

cat */*/*.introns.tsv | sed  -e 's/\.\./\t/' | grep GoldMap | grep -v GOLD | sort -k 2,2 -k 3,3 | gawk -F '\t' -f scripts/introns_precision_recall.awk  | sort -u | sed  -e 's/t1r/T1R/g'  -e 's/t2r/T2R/g'  -e 's/t3r/T3R/g' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}'  | sed -e 's/PFAL/Malaria/g' -e 's/STARlongzz/STAR long/g' >> $toto

set toto=RESULTS/Mapping_accuracy.light.txt
echo -n "## $toto :" > $toto
date >> $toto 
cat scripts/mapping_accuracy.header   >> $toto

cat */*/*.introns.tsv | grep -v R2 | grep -v R3 |  sed  -e 's/R1//g' | grep -v SRR | grep -v Simulated | sed  -e 's/\.\./\t/' | grep GoldMap | grep -v GOLD | sort -k 2,2 -k 3,3 | gawk -F '\t' -f scripts/introns_precision_recall.awk  | sort -u | sed  -e 's/r1//g'   | grep -v r2 | grep -v r3 | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto

# goto phaseLoop

##### INTRONS report  Insertion Deletion

foreach type (Intron  Insertion Deletion)

  set toto1="RESULTS/$type"_per_coverage_stats.txt
  set toto1L="RESULTS/$type"_per_coverage_stats.light.txt
  set toto3="RESULTS/$type"_per_coverage_stats.best.txt
  set toto4="RESULTS/$type"_per_coverage_stats.1support.txt
  set toto3L="RESULTS/$type"_per_coverage_stats.best.light.txt
  set toto4L="RESULTS/$type"_per_coverage_stats.1support.light.txt
  echo -n "### $toto1 :" > $toto1
  echo -n "### $toto3 :" > $toto3
  echo -n "### $toto4 :" > $toto4
  echo -n "### $toto1L :" > $toto1L
  echo -n "### $toto3L :" > $toto3L
  echo -n "### $toto4L :" > $toto4L
  date >> $toto1
  date >> $toto3
  date >> $toto4
  date >> $toto1L
  date >> $toto3L
  date >> $toto4L

 if ($type == Intron) then
    echo "## An alignment supporting an intron is defined by a line in the SAM/BAM file where the CIGAR contains an N with minimal accepted intron length 50 bases" > $toto1.1
    echo "## When a read is aligned at multiple sites, each of its alignments supporting an intron is counted" >> $toto1.1
    echo "## Some spliced genes are truly repeated, some are very similar. If one rejects all multiply aligned reads, the introns of these genes cannot be detected," >> $toto1.1
    echo "## Therefore, we keep the introns detected by multiply aligned reads, and do not artificially overestimate the specificity of methods unable to select the true positions" >> $toto1.1
    echo "## Note that in the benchmark, all reads are uniquely aligned, yet some support multiple neighboring introns" >> $toto1.1
 else
    echo "## An alignment supporting an indel is defined by a line in the SAM/BAM file where the CIGAR contains an I or o D" > $toto1.1
    echo "## When a read is aligned at multiple sites, each of its alignments supporting an indel is counted" >> $toto1.1
  endif

  echo -n "# Species\tRun\tMethod\tMinimal $type support" >> $toto1.1
  echo "\t$type in benchmark\t$type discovered in method\tFP: False positive $type\tTP: True positive $type\tFN: False negative $type\t$type discovery precision p=TP/(TP+FP)\t$type discovery recall r=TP/(TP+FN)\t$type discovery F-score 2pr/(p+r)" >> $toto1.1

  cat $toto1.1 >> $toto1
  cat $toto1.1 >> $toto1L
 
  if (-e $toto1.2) \rm $toto1.2

  foreach mm ($methods)
    cat $mm/*/*.delins.tsv | grep $type | sed -e 's/on_/on/' -e "s/$type//" | grep -v GOLD | sort -k 1,1n -k 2,2 |  sed  -e 's/\.\./\t/'  -e 's/g_//' | gawk '/^#/{next;}/GOLD/{next;}{support=$1;species="Human";run=$2;if(support > 200 && run != "Illumina")next;if(substr(run,1,4)=="PFAL")species="Malaria";method=$3;printf("%s\t%s\t%s\t%d",species,run,method,support); printf("\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", $5, $6,$8,$7,$9,$10,$11,$12);}'  >> $toto1.2
  end

  cat $toto1.2 | sort -k 1,1 -k 2,2 -k 3,3 -k 4,4n >  $toto1.3
  cat $toto1.3 | gawk -F '\t' '/^#/{print;next;}{if($3!=old)printf("\n\n");old=$3;print;}'  >> $toto1
  cat $toto1.3 | grep -v r2 | grep -v r3 | grep -v SRR | grep -v Simulated | gawk -F '\t' '/^#/{print;next;}{if($3!=old)printf("\n\n");old=$3;print;}' >> $toto1L

  echo "## Limited to best support-depth" >> $toto3
  cat $toto1.1 >> $toto3
   cat $toto1.3 | gawk -F '\t' 'BEGIN{best=0;}{z=$1 "\t" $2 "\t" $3; if (z != old) {old=z;if(bestScore)print best;bestScore=$12;best=$0;}if($12>bestScore){bestScore=$12;best=$0;}}END{if(bestScore)print best;}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto3
  
  echo "## Limited to best support-depth" >> $toto3L
  cat $toto1.1 >> $toto3L
   cat $toto1.3 | grep -v r2 | grep -v r3 | grep -v SRR | grep -v Simulated | gawk -F '\t' 'BEGIN{best=0;}{z=$1 "\t" $2 "\t" $3; if (z != old) {old=z;if(bestScore)print best;bestScore=$12;best=$0;}if($12>bestScore){bestScore=$12;best=$0;}}END{if(bestScore)print best;}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto3L
  
  echo "## Limited to 1 support" >> $toto4
  cat $toto1.1 >> $toto4
  cat $toto1.3 | gawk -F '\t' '{if($4==1)print}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto4

  echo "## Limited to 1 support" >> $toto4L
  cat $toto1.1 >> $toto4L
  cat $toto1.3 | grep -v r2 | grep -v r3 | grep -v SRR | grep -v Simulated | gawk -F '\t' '{if($4==1)print}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto4L


 #\rm $toto1.[123]
end
# goto phaseLoop  # fall through to aliLn

##############################################################################
##############################################################################
####### aliLn

phase_aliLn:

set toto=RESULTS/Aligned_length.histo.txt
echo -n  "## $toto : " > $toto
date >> $toto

\rm $toto.*
  foreach run (iRefSeq Roche PacBio Illumina SRR5189652 SRR5189667 SRR5437876 SRR5438850)

    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run.fastq.gz) then
        bin/dna2dna -i Fasta/$run/$run.fastq.gz -I fastq -getTM > Fasta/$run/$run.TM
      endif
    endif
    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run.fasta.gz) then
        bin/dna2dna -i Fasta/$run/$run.fasta.gz -I fasta -getTM > Fasta/$run/$run.TM
      endif
    endif
    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run'_1'.fasta.gz) then
        zcat  Fasta/$run/$run'_'?.fasta.gz | bin/dna2dna -I fasta -getTM > Fasta/$run/$run.TM
      endif
    endif

    set delta=1
    if ($run == Roche) set delta=10
    if ($run == PacBio) set delta=30
    if ($run == SRR5189652) set delta=30
    if ($run == SRR5189667) set delta=30
    if ($run == SRR5437876) set delta=30
    if ($run == SRR5438850) set delta=30
    if ($run == iRefSeq) set delta=100
    foreach mm ($methods)
      if (-e  $mm/$run/$mm.$run.aliLn) then
        cat $mm/$run/$mm.$run.aliLn | gawk '{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=$delta mm=$mm rr=$run >> $toto.1
      endif  
    end
  end 

  echo "Illumina\ttruth\t101\t217498656" >>  $toto.1
  echo "SRR5437876\ttruth\t300\t32935604" >>  $toto.1

  cat Fasta/iRefSeq/iRefSeq.TM | gawk -F '\t' '/^#/{next;}{k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=iRefSeq > $toto.1a
  cat Fasta/PacBio/PacBio.TM | gawk -F '\t' '/^#/{next;}{k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=30 mm="truth" rr=PacBio >> $toto.1a
  cat Fasta/SRR5189652/SRR5189652.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=30 mm="truth" rr=SRR5189652 >> $toto.1a
  cat Fasta/SRR5189667/SRR5189667.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=30 mm="truth" rr=SRR5189667 >> $toto.1a
  cat Fasta/Roche/Roche.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=10 mm="truth" rr=Roche >> $toto.1a



  echo "Human_T1\ttruth\t100\t20000000" > $toto.t0
  echo "Human_T2\ttruth\t100\t20000000" >>  $toto.t0
  echo "Human_T3\ttruth\t100\t20000000" >>  $toto.t0
  echo "Malaria_T1\ttruth\t100\t20000000" >>  $toto.t0
  echo "Malaria_T2\ttruth\t100\t20000000" >>  $toto.t0
  echo "Malaria_T3\ttruth\t100\t20000000" >>  $toto.t0
  
  set delta=1
  foreach tt (1 2 3)
    if (-e $toto.t$tt) \rm $toto.t$tt
    foreach mm ($methods)
      cat $mm/HG19t$tt'r1'/$mm.HG19t$tt'r1'.aliLn | gawk '{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=$delta mm=$mm rr=Human_T$tt >> $toto.t$tt
      cat $mm/PFALt$tt'r1'/$mm.PFALt$tt'r1'.aliLn | gawk '{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=$delta mm=$mm rr=Malaria_T$tt >> $toto.t$tt
    end

  end

  cat $toto.t[0123] >> $toto.1
  cat $toto.1 | sort -k 1,1 -k 3,3n -k 2,2 > $toto.2
  cat $toto.1a | sort -k 1,1 -k 3,3n -k 2,2 >> $toto.2

  echo -n "## $toto :" > $toto
  if (-e $toto.5) \rm $toto.5
  date >> $toto
  echo "## Histogram of length to be aligned in truth dataset,and aligned by each program. Each read is counted only once, at the location of its BAM primary alignment (excluding the secondaries with flag 256)" >> $toto
  foreach rr (iRefSeq Roche PacBio Illumina $pacbio_runs $long_illumina_runs)
    echo "# $rr\tTruth" > $toto.4
    cat $toto.2 | grep $rr | gawk -F '\t' '{k=$3;n=$4;m=$2;if(k>kMax)kMax=k;kk[k]=1;nn[m,k]=n;}END{kk[0]=1;for (k=0;k<=kMax;k++)if(kk[k]>0)printf("%d\t%d\n", k,nn["truth",k]);}' >> $toto.4
    echo "\n\n" >> $toto.4
    cat $toto.4 | scripts/transpose >> $toto.5
    echo "\n\n" >> $toto.5
  end
#
  foreach rr (iRefSeq Roche PacBio Illumina $pacbio_runs $long_illumina_runs)
    echo -n "# $rr\tTruth" > $toto.4
    foreach mm ($methods)
      echo -n "\t$mm" >> $toto.4
    end
    cat $toto.2 | grep $rr | gawk -F '\t' '{k=$3;n=$4;m=$2;if(k>kMax)kMax=k;kk[k]=1;nn[m,k]=n;}END{kk[0]=1;for (k=0;k<=kMax;k++)if(kk[k]>0)printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", k,nn["truth",k],nn["10_MagicBLAST",k],nn["20_HISAT2_relaxed",k],nn["21_HISAT2",k],nn["30_STAR",k],nn["31_STARlong",k],nn["32_STAR.2.6c",k],nn["40_TopHat2",k]);}' >> $toto.4
    echo "\n\n" >> $toto.4
    cat $toto.4 | scripts/transpose >> $toto.5
    echo "\n\n" >> $toto.5
  end

  foreach rr (Human_T1 Human_T2 Human_T3 Malaria_T1 Malaria_T2 Malaria_T3)
    echo -n "# $rr\tTruth" > $toto.4
    foreach mm ($methods)
      echo -n "\t$mm" >> $toto.4
    end
    cat $toto.2 | grep $rr | gawk -F '\t' '{k=$3;n=$4;m=$2;if(k>kMax)kMax=k;kk[k]=1;nn[m,k]=n;}END{kk[0]=1;for (k=0;k<=kMax;k++)if(kk[k]>0)printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", k,nn["truth",k],nn["10_MagicBLAST",k],nn["20_HISAT2_relaxed",k],nn["21_HISAT2",k],nn["30_STAR",k],nn["31_STARlong",k],nn["32_STAR.2.6c",k],nn["40_TopHat2",k]);}' >> $toto.4
    echo "\n\n" >> $toto.4
    cat $toto.4 | scripts/transpose >> $toto.5
    echo "\n\n" >> $toto.5
  end

  cat $toto.5 | scripts/transpose | sed -e 's/Truth/Actual reads/g' >> $toto

  cat RESULTS/Intron_per_coverage_stats.txt | gawk -F '\t' '/^#/{next}{z=$1 "\t" $2 "\t" $3;k=0+n[z];if($4>k)n[z]=$4;}END{for(k in n)printf("%s\t%d\n",k,n[k]);printf("toto\n");}' | sort | gawk -F '\t' 'BEGIN{printf("# Species\tRun\tMethods\tMaximal intron support\n")}{z=$1 "\t" $2;if (z != oldz){if(length(oldz)>3)printf("%s\t%s\t%s\n",oldz,substr(m[oldz],2),substr(n[oldz],2));m[z]="";n[z]="";}oldz=z;m[z]=m[z]","$3;n[z]=n[z]","$4;}' > RESULTS/Intron_per_coverage_stats.title.txt
\rm $toto.*

goto phaseLoop

##############################################################################
##############################################################################
## Done

phase6:
phaseLoop:
  echo -n "$1 done : "
  date

##############################################################################
##############################################################################
##############################################################################
##############################################################################
 
