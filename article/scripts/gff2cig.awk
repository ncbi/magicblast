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
#

{ 
    if ($3 != "exon") 
	next ;
}
{
    split ($9,aa,"Genbank:") ; 
    split(aa[2],bb,",") ; split(bb[1],cc,";") ; seq=cc[1] ;
    
    if(substr(seq,1,2) != "NM" && substr(seq,1,2) != "zNR") 
	next ;

    seq=seq ":" $1 ; chrom[seq] = $1 ; 
    nx[seq]++ ; i=nx[seq] ;
    a1[seq,i] = $4 ; a2[seq,i] = $5 ; strand[seq]=$7;
}
END {
    for (seq in nx)
    {
	n = nx[seq] ;
	printf ("%s\t%s",seq,chrom[seq]) ;
	if (strand[seq] == "+")
	{
	    printf ("\t%d\t%d\t", a1[seq,1], a2[seq,n]) ;
	    for(i = 1 ; i <=n ; i++)
	    {
		if (i>1) 
		{
		    dx = a1[seq,i] -a2[seq,i-1] - 1 ;
		    printf("%dN",dx);
		}
		dx = a2[seq,i] - a1[seq,i] + 1;
		printf ("%dM",dx);
	    }
	}
	else 
	{
	    printf ("\t%d\t%d\t", a1[seq,n], a2[seq, 1]) ;
	    for (i = n ; i >=1 ; i--)
	    {
		if (i<n)
		{
		    dx = a1[seq,i] - a2[seq,i+1] - 1 ;
		    printf("%dN",dx);
		}
		dx = a2[seq,i] - a1[seq,i] + 1;
		printf("%dM",dx);
	    }
	}
	printf("\t.\t%s\t.\n",strand[seq]);
    }
}



