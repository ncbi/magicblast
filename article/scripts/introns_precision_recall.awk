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

/^#/ {
    tt++ ;
    if(tt==1)
	print ;
    next ; 
}
{ 
    z1 = substr($2,1,4) ;
    z2 = substr($2,5) ;
    printf ("%s\t%s\t%s", z1, z2, $3) ; 
    
    u = $14 ; nu = $15 ;
    printf ("\t%s\t%d\t%d\t%d", $5, u+nu,u,nu) ; 
    
    c = 0 ; if (nu>0) c = nu / (u+nu) ;
    printf ("\t%.4f", c) ; 

    fp1 =$12 ; fp2 = $13 ; fp = fp1 + fp2 ;
    tp1 =$10 ; tp2 = $11 ; tp = tp1 + tp2 ;
    fn = $16 ;
    printf ("\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d", fp,fp1,fp2, tp1,tp2, tp, fn) ; 
    
    p = 0 ; r = 0 ; f = 0 ; c = 0 ;
    if (tp > 0)
    {
	p = tp / (tp + fp) ;
	r = tp / (tp + fn) ;
	f = 2 * p * r / (p+r) ;
	c = tp2 / tp ;
    }
    printf ("\t%.4f\t%.4f\t%.4f\t%.4f",p,r,f,c) ; 

    fp1 = $8 ; fp2 = $9 ; fp = fp1 + fp2 ; 
    tp1 = $6 ; tp2 = $7 ; tp = tp1 + tp2 ;
    fn=$16; 
    c = 0 ; if (tp > 0)	c = tp1 / tp ;

    printf ("\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d", fp, fp1, fp2, tp1, tp2, tp, fn) ;
    
    p = 0 ; r = 0 ; f = 0 ; c = 0 ;
    if (tp > 0)
    {
	p = tp / (tp+fp); 
	r = tp / (tp+fn);
	f = 2 * p * r / (p+r) ;
	c = tp2 /tp ;
    }
    printf ("\t%.4f\t%.4f\t%.4f\t%.4f",p,r,f,c) ; 
    
    printf("\n");
}
