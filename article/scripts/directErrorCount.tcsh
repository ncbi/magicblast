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
#  Author: Jean Thierry-Mieg  mieg@ncbi.nlm.nih.gov
#
## Direct statistics of the error counts reported in the BAM files
## We rely on the presence in the BAM files of the NM:i:x
## optional filed, collate the X values and report the statistics

set mm=$1
set run=$2

  set nreads=`cat Fasta/$run/*.count  | gawk '/^Fragment_kept/{n+=$2}END{print n}'`

zcat $mm/$run/$mm.$run.sam_sorted.gz | gawk -F '\t' -f scripts/directErrorCount.awk | scripts/tags | gawk '{split ( $1,aa,":" ) ; x=aa[3] ; y = $2 ; if(x>=0) ali+= y ; printf ( "%d\t%d\n",x,y ) ; }END{printf ("-5\t%d\n",nreads - ali)}'  nreads=$nreads | sort -k 1,1nr | gawk '{if (line<1)print "-6:"method": -5:unaligned -4:partial with error,  -3:partial no error, -2: complete with error, -1: complete no errors, 0: no error, 1,2,3...:n errors, columns 3 and 4 are cumuls";line++;if($1>=0){n1 += $2 ; n2 += $1*$2 ;}printf ( "%d\t%d\t%d\t%d\n",$1,$2,n1,n2 ) ; }' method=$mm nreads=$nreads | sort -k 1,1n >  $mm/$run/$mm.$run.nerrors

