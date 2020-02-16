#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mcheck.h>
#include <time.h>
/* 
#include <gmp.h> 
#include <mpfr.h>
#include <mpc.h>
*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_cdf.h>
#include </home/BJARKI/verk/truncation/wfiles/statistics.h>
/* Motivation for Kolmogorov distance between distributions */
/* Kolmogorov distance, Kantorovich-Rubinstein distance, the bounded-Lipschitz distance, Cramer-von Mises distance, Wassertein metric, Anderson-Darling */
static const double CUTOFF_K = 0.;
static const double ALPHA = 1.;
static const double K = 0.;
static const double NUMBRUNS = 10.;

enum{
  SAMPLESIZE = 158};
enum{
  LINKAGEGROUPS = 20};



gsl_rng *rngtype;


static void setup_rng(unsigned long int s)
{
const gsl_rng_type *T;
gsl_rng_env_setup();
T= gsl_rng_default;
rngtype= gsl_rng_alloc(T);
gsl_rng_set(rngtype,s);
}


static void printmatrix( gsl_matrix *M, int radir, int dalkar )
{
  int i,j;
  for( i = 1; i <= radir; i++){
    for( j = 1; j <= dalkar; j++){
      printf("%g ", gsl_matrix_get(M, i,j));}
    printf("\n");}
}



static void initmatrix( gsl_matrix_int *G )
{
  /* n is number of chroms per linkage group; h is number of linkage groups */
  int i,j;
  for( j = 1; j <= LINKAGEGROUPS ; j++){
    gsl_matrix_int_set(G,j, 0,  SAMPLESIZE);
    assert( gsl_matrix_int_get(G,j,0) == SAMPLESIZE);
    for( i = 1; i <= SAMPLESIZE ; i++){
      gsl_matrix_int_set(G, j,i, i);}}

 /* set the number of linkage groups  */
  gsl_matrix_int_set(G,0,0, LINKAGEGROUPS );
}


static int assigntogroup( double x, gsl_rng *r)
{
  int  group = 0;
  double u ;
  if( gsl_rng_uniform_pos(r) <= x ){
    /* given block participates in large family */
    u = gsl_rng_uniform(r);
    group = (u < .25 ? 1 : ( u < .5 ? 2 : (u < .75 ? 3 : 4)));
    assert( group > 0);
    assert( group < 5); }
  /* return the number of the parental chromosome; if chrom= 0 then block is not in large family */
  return( group );
}



static int checkindex(int i, int mgroup, gsl_matrix_int *I)
{
  /* I is 4 \times n matrix of block indexes merging in each group */
  /* I collects indexes of blocks merging into merging groups for each linkage group */
  /* check if block with index i among blocks merging in group  mgroup*/
  int j = 1;
  int stop = 0;
  assert( i > 0);
 /* I[mgroup,0] is number of blocks in merging in merging group mgroup */
  assert( gsl_matrix_int_get(I, mgroup, 0) > 0);
    /* stop = ((i == gsl_matrix_int_get(I, mgroup, 1)) || (i == gsl_matrix_int_get(I, mgroup, gsl_matrix_int_get( I, mgroup,0)))) ? 1 : 0; */
  if( i <= gsl_matrix_int_get(I, mgroup, gsl_matrix_int_get( I, mgroup,0))){
      do{
	stop = (gsl_matrix_int_get( I, mgroup,j) == i ? 1 : 0);
	j = j + 1;}
      while( (j <= gsl_matrix_int_get(I, mgroup,0) ) && (stop == 0));}
  
  assert( (stop == 0) || (stop == 1));
  return( stop);
}


static int samplepairwise( gsl_matrix_int *G,  double *pairlocus, gsl_rng * r)
{
  /* sample locus for a  pairwise merger */
  int j;
  for(j = 1; j <= LINKAGEGROUPS; j++){
    pairlocus[j] = (gsl_matrix_int_get(G, j, 0) > 1 ? gsl_sf_choose( gsl_matrix_int_get(G, j, 0), 2) : 0.);}
  gsl_ran_discrete_t * P;
  P = gsl_ran_discrete_preproc(LINKAGEGROUPS + 1, pairlocus);
  int svar = (int)gsl_ran_discrete( r, P);
  gsl_ran_discrete_free(P);
  assert( gsl_matrix_int_get(G, svar, 0) > 1 );
  assert( svar >= 1);
  assert( svar <= LINKAGEGROUPS);
  /* return the linkage group at which two blocks merge */
 return( svar);
}


static void picktwoblocks(gsl_matrix_int *G,  int lgroup, int *indexes,  int *blocks,   gsl_rng *r)
{
  assert( gsl_matrix_int_get(G, lgroup, 0) > 1);
  int k = 1;
  int i = 0;

  blocks[0] = blocks[1] = 0;
  indexes[0] = indexes[1] = 0;
  /* collect all indexes of active blocks at linkage group lgroup */
  for( k = 1; k <=  SAMPLESIZE; k++){
    if( gsl_matrix_int_get(G, lgroup, k) == k){
    indexes[i] = k;
    i = i + 1;}}

  assert( i > 1);
  assert( indexes[0] > 0);
  assert( indexes[1] > 0);

  /* sample two block indexes from  indexes of active blocks at lgroup */
  gsl_ran_choose(r, blocks, 2, indexes, gsl_matrix_int_get(G, lgroup, 0), sizeof(int));
  assert(blocks[0] > 0);
  assert(blocks[1] > 0);
}


static double randombeta( gsl_rng *r)
{
  /*return a random variate from a complete beta */
  /* or an incomplete beta */
  /*  the  expected value m as N to infty is taken as  1./(\ALPHA - 1) */
  /* if K > 0, draw U uniform(0,  K/(m + K)); */
  /* then  X = F^{-1}(U)  where F is inverse Beta(2-ALPHA, ALPHA) */
  
   /* double u = gsl_ran_flat( r, 0.,  K/(K + 1./(ALPHA - 1.)));  */
  return( K > (0.) ? gsl_cdf_beta_Pinv( gsl_ran_flat( r, (0.),  K/(K + (1.)/(ALPHA - 1.))), (2.) - ALPHA, ALPHA) : gsl_ran_beta(r,  (2.) - ALPHA, ALPHA));
}



static int mergegroups(gsl_matrix_int *G, gsl_matrix_int *I, int *indexes,  double *plocus, int *pairblocks, double x, gsl_rng *r)
{
  /* lgroup is linkage group for which to merge blocks */
  /* double-check if still segregating blocks at lgroup */
  int k, i, j, g, merger, pairwiselocus; 
  merger = 0;
  /*  double x = ( ALPHA < 2. ? randombeta(r) : 0. );  */

  /* draw the linkage group seeing a pairwise merger */
  pairwiselocus = samplepairwise( G, plocus, r);

  /* draw two blocks at lgroup picked for pairwise merger */
  picktwoblocks(G, pairwiselocus, indexes,   pairblocks, r);

  for( j = 1; j <= LINKAGEGROUPS; j++){
    if( gsl_matrix_int_get( G, j, 0) > 1){
  gsl_matrix_int_set( I, 1,0, 0);
  gsl_matrix_int_set( I, 2,0, 0);
  gsl_matrix_int_set( I, 3,0, 0);
  gsl_matrix_int_set( I, 4,0, 0);
  /* mark the locus with the pairwise merger 
  gsl_matrix_int_set( I, 1, 0, (pairwiselocus == j ? 2 : 0)); */
  /**/

  for( k = 1; k <= SAMPLESIZE; k++){
    if( gsl_matrix_int_get( G, j, k) == k){
      /* active block with index k; assign to group */
      /* g = ( (j == pairwiselocus) && ((k == pairblocks[0]) || (k == pairblocks[1])) ? assigntogroup(1., r) : assigntogroup(x, r)); */
      /* assign the two blocks picked initially to merger group 1; assures a merger */
       g = ( (j == pairwiselocus) && ((k == pairblocks[0]) || (k == pairblocks[1])) ? 1 : assigntogroup(x, r));
      /* update number of blocks assigned to same merging group */
    gsl_matrix_int_set(I, g, 0, gsl_matrix_int_get(I, g, 0) + (g > 0 ? 1 : 0));
    gsl_matrix_int_set(I, g, gsl_matrix_int_get(I, g, 0), (g > 0 ? k : gsl_matrix_int_get(I, g, gsl_matrix_int_get(I, g,0))));}}
     /* participates in family and  assigned to group */
     /* check if a merger occurs at linkage group j */
  merger = merger +  (gsl_matrix_int_get(I,1,0) > 1 ? 1 : 0) + (gsl_matrix_int_get(I,2,0) > 1 ? 1 : 0) + (gsl_matrix_int_get(I,3,0) > 1 ? 1 : 0) + (gsl_matrix_int_get(I,4,0) > 1 ? 1 : 0);

  /* lines = lines +  (gsl_matrix_int_get(I,1,0) > 1 ? gsl_matrix_int_get(I,1,0) : 0) + (gsl_matrix_int_get(I,2,0) > 1 ? gsl_matrix_int_get(I,2,0)  : 0) + (gsl_matrix_int_get(I,3,0) > 1 ? gsl_matrix_int_get(I,3,0) : 0) + (gsl_matrix_int_get(I,4,0) > 1 ? gsl_matrix_int_get(I,4,0) : 0);
   */
    for( i = 1 ; i <= 4; i++){
      if( gsl_matrix_int_get(I, i, 0) > 1){
	assert( gsl_matrix_int_get(I, i, 1) < gsl_matrix_int_get(I, i, 2));
	assert( gsl_matrix_int_get(I, i, 2) <= gsl_matrix_int_get(I, i, gsl_matrix_int_get(I, i, 0)));
	assert( gsl_matrix_int_get(I, i, 1) < SAMPLESIZE);
	/* at least two lines merge in group i */
	for( k = gsl_matrix_int_get(I, i, 1); k <= SAMPLESIZE ; k++){
	  assert( k > 0);
	  assert( j > 0);
	  gsl_matrix_int_set(G, j, k, checkindex( gsl_matrix_int_get(G, j, k), i, I) == 1 ? gsl_matrix_int_get(I, i, 1) : gsl_matrix_int_get(G,j, k));}
    /* update the count of active blocks at linkage group j */
	gsl_matrix_int_set(G, j, 0, gsl_matrix_int_get(G,j,0) - (gsl_matrix_int_get(I,i,0) > 1 ? (gsl_matrix_int_get(I,i,0) - 1) : 0));
	assert( gsl_matrix_int_get(G, j, 0) > 0);
	assert( gsl_matrix_int_get(G, j, 0) <= SAMPLESIZE);
      }}
    /* update number of  linkage groups with active  loci */
    gsl_matrix_int_set(G, 0, 0, gsl_matrix_int_get(G,0,0) - (gsl_matrix_int_get(G,j,0) == 1 ? 1 : 0));} }
  /* return if merger > 0 then at least one pairwise merger overall */
  /* run until at least one pairwise merger overall */
  
  /* assert( (merger > 0 ? lines == 2 : lines == 0) ); */
  return( merger);
}


static double numerator( double x, gsl_matrix_int *G)
{
   int j;
   double y = 1.;
   double n;
   /* check at least one locus still with blocks */
   assert( gsl_matrix_int_get(G, 0, 0) > 0 );
   for( j = 1; j <= LINKAGEGROUPS; j++){
     n = (double)gsl_matrix_int_get(G, j,0);
     /*  y = y*( (n > 1. ? (pow( 1. - x, n) +  (n*x*pow(1.-x, n-1.))) : 1.));} */
       
     y =  y * (gsl_matrix_int_get(G,j,0) > 1 ? ( (pow(1-x,n)) + (4.*n*(x/4.)*pow(1.-x, n-1.)) +  (4.*3.*gsl_sf_choose( gsl_matrix_int_get(G,j,0),2)*pow(x/4.,2.)*pow(1.-x,n-2.)) + (gsl_matrix_int_get(G,j,0) > 2 ? (4.*3.*2.*gsl_sf_choose( gsl_matrix_int_get(G,j,0),3)*pow(x/4,3)*pow(1.-x,n-3.)) : 0.) +  (gsl_matrix_int_get(G, j,0) > 3 ? (4.*3.*2.*gsl_sf_choose(n,4.)*pow(x/4.,4.)*pow(1.-x,n-4.)) : 0.)) : 1.); }

     /* return the f function */
   return( 4.*(1. - y)/(x*x) );
}


static double sampletime( gsl_matrix_int *G,  double *y,  gsl_rng *r )
{
  int j;
  double rate, ratee,   nm, exactf;
  double timi = 0. ;
  rate = 0.;
  ratee = 0.;
  for( j = 1; j <= LINKAGEGROUPS; j++){
    rate = rate  +  (gsl_matrix_int_get( G, j, 0) > 1 ? gsl_sf_choose( gsl_matrix_int_get( G, j, 0),2) : 0);
    ratee = ratee + (gsl_matrix_int_get( G, j, 0) > 2 ? gsl_sf_choose( gsl_matrix_int_get( G, j, 0),3) : 0);}
  assert( rate > 0. );
  
  do{
    timi = timi +  gsl_ran_exponential(r, 1./rate);
    y[0] = randombeta(  r );
    exactf = numerator( y[0], G) ;
    nm = (y[0] < 0.0001 ? (rate  - (2.*y[0]*ratee)) : exactf)/rate ;
    /* printf( "%g %g\n", x, nm); */
    /* assert( nm <= 1); */
  }
  while( gsl_rng_uniform(r) > (ALPHA < 2. ? (nm) : .25) );
  assert( timi > 0. );
  return( timi);
}



static void genealogy(  gsl_matrix_int *G, gsl_matrix_int *I,  int *indexes, double *probl,   int *bl,  gsl_matrix *BLS, double *y,  gsl_rng *r)
{
  /* run a genealogy */
  /* n is sample size, ie number of chromosomes per  linkage groups */
  /* h is number of linkage groups */
  /* G[0,0] is the number of groups with  segregating blocks, ie not reached mrca */
  /* I is a 4 \times n matrix, each row holds the indexes of blocks merging in each group */

  /*x is prob of being in big family */
  
  initmatrix(G);
 double timi = 0.;
 int merger = 0; 
 /* clock_t begin, end; */

  gsl_matrix_set_zero(BLS);
  /*
  timi =  sampletime( G, r);
  updatebls(SAMPLESIZE, LINKAGEGROUPS, G, BLS, timi);
  */
  while( gsl_matrix_int_get(G,0,0) > 0){
    /* at least 1 linkage group with segregating blocks */ 
    /* draw time and merge groups */
    timi =  sampletime(G, y,  r);
    updatebls(SAMPLESIZE, LINKAGEGROUPS, G, BLS, timi);
    /* update forest */
    merger =  0;
    do{
      /* begin = clock(); */
      merger =  mergegroups(G, I, indexes,  probl, bl, y[0],  r);
      /* end = clock();
	 printf("%g\n", ((double)(end - begin)) / ((double)CLOCKS_PER_SEC) ); */
     }
    while( merger < 1); }
  /* obtained a new configuration of forest */
}




static void runstatistics( long rseed)
{
  /*run many genomealogies and record statistics */
  double b = 0.;
  
   gsl_matrix_int *G = gsl_matrix_int_calloc( LINKAGEGROUPS + 1,  SAMPLESIZE + 1);
   gsl_matrix_int *I = gsl_matrix_int_calloc( 5,  SAMPLESIZE + 1);
   gsl_matrix *BLS = gsl_matrix_calloc( LINKAGEGROUPS + 1,  SAMPLESIZE + 1);
   gsl_matrix *CBLS = gsl_matrix_calloc( LINKAGEGROUPS + 1,  SAMPLESIZE + 1);
   double *problocus = (double *)calloc( LINKAGEGROUPS + 1, sizeof(double));
   int *pairblocks = (int *)calloc( 2, sizeof(int));
  int *indexes = (int *)calloc( SAMPLESIZE, sizeof(int));
  double * coords1 = (double *)calloc( LINKAGEGROUPS + 1, sizeof( double));
  double * coords2 = (double *)calloc( LINKAGEGROUPS + 2, sizeof( double));
    double *y = (double *)malloc(sizeof(double));
    clock_t begin, end;
   setup_rng( rseed );
  while(b < NUMBRUNS){
    /**/
    // begin = clock();
    genealogy(G, I,  indexes,   problocus, pairblocks, BLS, y, rngtype);
    /* printmatrix( BLS, LINKAGEGROUPS, SAMPLESIZE ); */
    /* printf("%g\n", gsl_matrix_get(BLS, 1,2) + gsl_matrix_get(BLS, 1,3) + gsl_matrix_get(BLS, 1,4) ); */
    // end = clock();
    /* printf("%g\n", ((double)(end - begin)) / ((double)CLOCKS_PER_SEC) ); */
    /* sumbls(SAMPLESIZE, LINKAGEGROUPS, BLS, coords); */
    /*static void collectbls(int n, int h, gsl_matrix *BLS, gsl_matrix *collectedbls )*/
    /* collectbls(SAMPLESIZE, LINKAGEGROUPS, BLS, CBLS); */
    /* updatesbls( SAMPLESIZE, LINKAGEGROUPS, G, SBS ); */
    /* printf( "%g %g\n", coords[0], coords[1]); */
    /* static void sfsarray(int n, int H,  gsl_matrix * BLS, double *coord1, double *coord2) */
    
     sfsarray(SAMPLESIZE, LINKAGEGROUPS, BLS, coords1, coords2);
     coords1[0] = 0.;
     coords2[0] = 0.;
    
     /* shift the data for correlation 
     for( j = 1; j <=  LINKAGEGROUPS ; j++){
       coords1[j-1] = coords1[j];
       coords2[j-1] = coords2[j]; } */
     /* print the corelations 
	printf("%g\n",  gsl_stats_correlation( coords1, 1, coords2, 1,  LINKAGEGROUPS)); */
     
     /*  print average of  coordinates over loci  */
     int j;
     for( j = 1; j <=  LINKAGEGROUPS ; j++){
       coords1[0] = coords1[0]  +  coords1[j] ;
       coords2[0] = coords2[0]  +  coords2[j] ;}
     printf("%g %g\n", coords1[0]/((double)LINKAGEGROUPS),  coords2[0]/((double)LINKAGEGROUPS));
     
     b = b + (1.0f);}

  /* static void printm( int n, int h, double denom,  gsl_matrix *X) */
  /*printm( SAMPLESIZE, LINKAGEGROUPS, NUMBRUNS,  CBLS); */
  
  /* free mem */
   gsl_matrix_int_free( G);
   gsl_matrix_int_free( I);
   gsl_matrix_free(BLS);
   gsl_matrix_free(CBLS);
   
   free( problocus);
   free( pairblocks);
   free(indexes);
   free(coords1);
   free(coords2);
   free(y);
}





 int main(int argc, char *argv[])
 { 
  

    runstatistics( atol(argv[argc-1]));
  
   
gsl_rng_free(rngtype);
  
   
   return GSL_SUCCESS;
 }
