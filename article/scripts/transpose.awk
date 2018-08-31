{ j++ ; for (i = 1 ; i <= NF ; i++)a[i,j] = $i ; if(NF>iMax) iMax = NF ; jMax = j ;}
END { 
  for (i = 1 ; i <= iMax ; i++) 
    {
      printf ("%s", a[i,1]) ;
      for (j = 2 ; j <= jMax ; j++)  
	printf ("\t%s", a[i,j]) ;
      printf ("\n",i) ;
    }
}

