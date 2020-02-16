/* include with genomic.c */



static unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}



static void printm( int n, int h, double denom,  gsl_matrix *X)
{

  int i, j;
  for( j = 1; j <= h; j++){
    for( i = 0; i < n; i++){
      printf("%g ", gsl_matrix_get(X, j, i)/denom);}
    printf( "\n");}
}





static void updatebls( int n, int h, gsl_matrix_int *G, gsl_matrix *B, double timi)
{
  /* update the branch length spectrum */
  /*h is linkage groups, n is sample size  */
  int i, j, k, size, blocks;
  for(j = 1; j <= h; j++){
    if ( gsl_matrix_int_get(G, j, 0) > 1){
      /* at least two blocks at lgroup j */
      blocks = 0;
    for( i = 1 ; i <= n; i++){
      if( gsl_matrix_int_get(G, j, i) == i){
	/* i is an index of an active block */
	size = 0;
	/* add to count of active blocks at lgroup j */
	blocks = blocks + 1;
	/* measure the size of the active block with index i*/
	for( k = i ; k <= n; k++){
	  size = size + (gsl_matrix_int_get(G, j, k) == i ? 1 : 0);}
	assert( size > 0);
	assert( size < n);
        /* add timi to branch length of size only if still active blocks at lgroup j */
	gsl_matrix_set(B, j, size, gsl_matrix_get(B, j,size) +  (gsl_matrix_int_get(G, j,0) > 1 ? timi : 0.));
	/* update total tree size excluding  external branches and  n-1 size */
	gsl_matrix_set(B, j, 0,  gsl_matrix_get(B,j,0) +  (gsl_matrix_int_get(G,j,0) > 1 ? (size > 1 ? (size < n-1 ? timi : 0.) : 0.) : 0.));}}
	/* update total tree size including  external and n-1 sizes 
	   gsl_matrix_set(B, j, 0,  gsl_matrix_get(B,j,0) + (gsl_matrix_int_get(G,j,0) > 1 ? timi : 0.));}}  */
    /* check count of active blocks */
    assert( blocks == gsl_matrix_int_get(G,j,0));
    }
  }
    /* update the total tree size at linkage group j; G[j,0] is number of blocks at lgroup j
       gsl_matrix_set(B, j, 0, gsl_matrix_get(B,j,0) +  (gsl_matrix_int_get(G,j,0) > 1 ? (((double)gsl_matrix_int_get(G, j,0)) * timi) : 0)); */
    /**/
}










static void sumbls(int n, int h, gsl_matrix *BLS, double *s)
{
  int j, i;
  double tail = 0.;
  s[0] = s[1] = 0. ;
  for( j = 1; j <= h ; j++){
    /* BLS[j,0] is the total tree size  of tree at linkage group j */
    s[0] = s[0] +  (gsl_matrix_get(BLS, j, 2) + gsl_matrix_get(BLS, j, 3) + gsl_matrix_get(BLS, j, 4))/gsl_matrix_get(BLS,j,0);
    tail = 0. ;
    for( i = 15; i <= n; i++){
      /*sum the tail for the SFS at linkage group j*/
      tail = tail +  (gsl_matrix_get(BLS, j, i)/gsl_matrix_get(BLS, j, 0));}
    s[1] = s[1] + tail;}
}


static void sfsarray(int n, int H,  gsl_matrix * BLS, double *coord1, double *coord2)
{
  /* sum over only certain sizes of the BLS */
  /* sums the relative branch lengths for each linkage group */
  /* check about differentials  D_i :=  R_i - R_{i+1} */ 
  int j, i;
  for( j = 1; j <= H ; j++){
    /* assert( gsl_matrix_get( BLS, j, 1) <= gsl_matrix_get( BLS,j, 0)); */
    coord1[j] = 0.;
    /* sum only 2nd, 3rd, 4th size classes for the 1st coordinate of the 2d statistic */
    /* for( i = 2; i <= 4; i++){ */
       coord1[j] = (gsl_matrix_get(BLS, j, 2) + gsl_matrix_get(BLS, j, 3) +  gsl_matrix_get(BLS, j, 4))/gsl_matrix_get(BLS,j,0); 
    /*  only relative   external branches  
	coord1[j] = (gsl_matrix_get(BLS, j, 1))/gsl_matrix_get(BLS,j,0); */
    coord2[j] = 0.;
    /* sum size classes  from  15 to  n-2 for the 2nd  coordinate of the 2d statistic */
    /* sum size classes from 15 to n-1 */
    for( i = 15; i < n-1; i++){
      assert( gsl_matrix_get( BLS, j, i) <= gsl_matrix_get( BLS,j, 0));
    coord2[j] = coord2[j] +   (gsl_matrix_get(BLS, j, i)/gsl_matrix_get(BLS, j, 0));}}
}




static void sfsarrayall(int n, int H,  gsl_matrix * BLS, double *coord1, double *coord2)
{
  /* sum over only certain sizes of the BLS */
  /* sums the relative branch lengths for each linkage group */
  /* check about differentials  D_i := 1(R_i > 0)*(R_i - R_{i+1})/R_i */ 
  /* include  both singletons and n-1 class */
  int j, i;
  for( j = 1; j <= H ; j++){
    coord1[j] = 0.;
       coord1[j] = (gsl_matrix_get(BLS, j, 1) + 0.)/gsl_matrix_get(BLS,j,0);
    /*  only relative   external branches  
	coord1[j] = (gsl_matrix_get(BLS, j, 1))/gsl_matrix_get(BLS,j,0); */
    
    coord2[j] = 0.;
    /* sum size classes  from  15 to  n-2 for the 2nd  coordinate of the 2d statistic */
    /* sum size classes from 15 to n-1 */
    for( i = 15; i < n ; i++){
    coord2[j] = coord2[j] +   (gsl_matrix_get(BLS, j, i)/gsl_matrix_get(BLS, j, 0));}}
}





static void collectbls(int n, int h, gsl_matrix *BLS, gsl_matrix *collectedbls )
{
  int i, j;
  for( j = 1; j <= h ; j++){
    /* BLS[j,0] is the total tree size  of tree at linkage group j */
    assert( gsl_matrix_get( BLS, j, 0) > 0. );
    for( i = 1; i < n ; i++){
      gsl_matrix_set( collectedbls, j, i, gsl_matrix_get(collectedbls, j, i)  +  (gsl_matrix_get( BLS, j, i)/gsl_matrix_get(BLS, j, 0)));}
    gsl_matrix_set( collectedbls, j, 0, gsl_matrix_get( collectedbls, j, 0) +  gsl_matrix_get(BLS, j, 0)); }
}


static void updatesbls( int n, int h, gsl_matrix_int *G, gsl_matrix *BS, double timi )
{ 

 /* update the sequence specific branch length spectrum */
  /*h is linkage groups, n is sample size  */
  /* BS is n\times h \times  n  matrix */
  int i, j, k, size;
  for(j = 1; j <= h ; j++){
    for( i = 1 ; i <= n; i++){
      if( gsl_matrix_int_get( G, j,i) == i){
	/* i is an index of an active block */
	size = 0;
	/* measure the size of the active block with index i*/
	/* size includes the index; size = 1 then external branch */
	/* size = 2 then block contains one other leaf besides i */
	for( k = i ; k <= n; k++){
	  size = size + (gsl_matrix_int_get(G, j, k) == i ? 1 : 0);}
	gsl_matrix_set(BS, (n*(j-1)) + i, size, gsl_matrix_get(BS, (n*(j-1)) + i, size) + timi); }}}
}


static void sumoverseqs(int n, int h, gsl_matrix *SBL, double *coord1, double *coord2)
{
  /* given a matrix of sequence specific  branch lengths, sum over  */
  int i, j, k;
  double tail = 0.;
  for( j = 1 ; j <= h; j++){
    coord1[j] = coord2[j] = 0.;
    /* run over the n sequences for linkage group j */
    for( i = 1; i <= n; i++){
      /* coord 1 is the sum over size classes  2, 3, 4 ; */
      coord1[j] = coord1[j] +  gsl_matrix_get( SBL, (n*(j-1)) + i, 2) + gsl_matrix_get( SBL, (n*(j-1)) + i, 3) + gsl_matrix_get( SBL, (n*(j-1)) + i, 4);
      tail = 0.;
      /* sum the tail minus last coordinate  for  sequence i */
      for(k = 15; k < n ; k++){
	tail = tail  +  gsl_matrix_get( SBL,  (n*(j-1)) + i, k);}
      /* add the tail to the sum over leaves */
      coord2[j] = coord2[j] + tail; } }
}
